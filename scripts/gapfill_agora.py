import pandas as pd
from framed import *
from framed.model.environment import *
from carveme.reconstruction.gapfilling import gapFill
import glob
import re
from framed.cobra.reconstruction import extend_model_with_DB_plainText

filter_model_id = "Clostridium_perfringens_ATCC_13124.xml"

output_path = "../report/tables/S8_gapfills_new.tsv"
studied_mediums = ["M3_M4", "M7", "M5", "M11", "M2", "M1", "M13", "M14"]
organisms = pd.read_table("../data/organisms.tab", comment="#")

# read growth matrix
growth = pd.read_table("../report/tables/S4_growth_matrix.tab")
growth["M3_M4"] = growth[["M3", "M4"]].max(axis=1)
growth = pd.melt(growth, id_vars=['Species'], var_name="Media", value_name="Viable")
growth.Viable = growth.Viable>0
growth = growth[growth.Media.isin(studied_mediums)]
growth = growth.merge(organisms[["agora", "species"]], left_on="Species", right_on='species')
viable_species = growth.query("agora == agora").groupby("agora").Viable.any()
viable_species = set(viable_species.index[viable_species])
print "Growth matrix loaded"

# Read medias
media_db = pd.read_table("../data/media_db.tab", comment="#")
media_db["compound"] = [cpd.strip() for cpd in media_db["compound"]]
media_db = media_db[["medium", "compound"]].groupby("medium").agg(lambda x: list(x)).to_dict()['compound']
mediums = {}
for m in studied_mediums:
    env = Environment.from_compounds(media_db[m], max_uptake=1000, exchange_format="'R_EX_{}_LPAREN_e_RPAREN_'")
    mediums[m] = env
print "Media database loaded"

# Read models
models = {}
simulations = {}
all_model_paths = glob.glob('../data/agora_f/*.xml')
for model_path in all_model_paths:
    model_name = os.path.basename(model_path)
    if model_name not in set(growth.agora.values):
        continue

    if model_name != filter_model_id and filter_model_id:
        continue

    print "{}/{} : Reading model '{}'".format(all_model_paths.index(model_path) + 1, len(all_model_paths), model_name)

    # with nostdout():
    model = load_cbmodel(model_path, flavor="fbc2")
    model.detect_biomass_reaction()

    biomass = model.reactions[model.biomass_reaction]
    if "M_spmd_c" in biomass.stoichiometry:
        del biomass.stoichiometry["M_spmd_c"]
    if "M_spmd_c" in biomass.stoichiometry:
        del biomass.stoichiometry["M_spmd_c"]

    models[model_name] = model
print "Read models"

# Load AGORA universe
universe = CBModel("agora_universe")
universe, db_reactions = extend_model_with_DB_plainText(universe, "../data/agora_reactions.bioopt")
print "Read universe model"

agora_reactions = pd.read_table("../data/agora_reactions.tab")
agora_reactions.EC[agora_reactions.EC.isnull()] = ""
agora_reactions.Subsystem[agora_reactions.Subsystem.isnull()] = ""
agora_reactions.ReactionName[agora_reactions.ReactionName.isnull()] = ""
print "Read universe annotation"


models_f = models #{m_id: model for m_id, model in models.iteritems() if m_id in viable_species}


biomass_target_flux = 0.1
models_gf = {}
for i, m_tuple in enumerate(models_f.iteritems(), start=1):
    m_id, model = m_tuple
    if m_id in viable_species:
        models_gf[m_id] = model.copy()
        print "{}/{}: Cloned model '{}'".format(i, len(models_f), m_id)
    else:
        models_gf[m_id] = model
        print "{}/{}: Skip cloning for model '{}'".format(i, len(models_f), m_id)
print "Cloned all models"

models_gf_sorted = list(sorted(models_gf.keys()))
print "Sorted all models"

#
# Start with writing header file and all reactions already existing in the model
#
with open(output_path, "w") as output_f:
    output_f.write("ModelID\tMediumID\tMediumNo\tReactionID\tEquation\tReactionName\tEC\tSubsystem\tGrowing\n")

#
# Loop through all media/species growth data combinations
#
for model_no, model_id in enumerate(models_gf_sorted, start=1):
    model = models_gf[model_id]  # .copy()

    print "{}/{} : {}".format(models_gf_sorted.index(model_id), len(models_gf_sorted), model_id)
    if model_id != filter_model_id and filter_model_id:
        continue

    #
    # Write reactions which were initially present in the model into output table
    #
    with open(output_path, "a") as output_f:
        for r_id, r in sorted(model.reactions.iteritems(), key=lambda x: x[0]):
            if re.search("biomass|growth", r_id, re.IGNORECASE):
                continue

            reaction_annotation = agora_reactions[agora_reactions.ReactionID == r_id]
            reaction_annotation = "\t".join(reaction_annotation[["ReactionName", "EC", "Subsystem"]].values[0])

            line = "\t".join([model_id, "", "", r.id, r.to_equation_string()])
            output_f.write("{}\t{}\t\n".format(line, reaction_annotation))

    for medium_no, medium_id in enumerate(studied_mediums):
        medium = mediums[medium_id].copy()
        medium.apply(models_gf[model_id], exclusive=True, warning=False)
        model = models_gf[model_id]

        # Get experimental growth
        is_growing_invitro = growth[(growth.Media == medium_id) & (growth.agora == model_id)].Viable.values[0]

        # Get simulation growth
        s = FBA(model)
        is_growing_insilico = s.message == "optimal" and s.fobj > biomass_target_flux or s.message == "unbounded"

        print "{}: {} {}".format(medium_id, is_growing_invitro, is_growing_insilico)
        #
        # Gapfill models where we know they should grow
        #
        reconstructed_model = None
        if is_growing_invitro and not is_growing_insilico:
            print "{:>5}\t{}".format(medium_id, model_id)
            # Call carveme and extract gapfilled reactions
            fixed_reactions = model.reactions.keys()
            biomass_eq = model.reactions[model.biomass_reaction].to_equation_string()

            reconstructed_model = gapFill(model, universe, constraints=medium, min_growth=biomass_target_flux,
                                          inplace=False, bigM=1e3)

            if not reconstructed_model:
                print "Gapfilling failed"
                with open(output_path, "a") as output_f:
                    line = "\t".join([model_id, medium_id, str(medium_no), r.id, r.to_equation_string()])
                    output_f.write("{}\t\t\t\t0\n".format(line, reaction_annotation))
                continue

            gapfilled_reactions = set(reconstructed_model.reactions.keys()) - set(model.reactions.keys()) - set(
                ["Growth"])
            gapfilled_reactions = {r.id: r for r in reconstructed_model.reactions.itervalues() if
                                   r.id in gapfilled_reactions}
            gapfilled_reactions_ids = gapfilled_reactions.keys()
            gapfilled_reactions_len = len(gapfilled_reactions)

            # Double check that all gap-filled reactions are really needed
            removed_reactions = []
            for r_id in gapfilled_reactions_ids:
                sol = reaction_deletion(reconstructed_model, [r_id], constraints={k: v for k, v in medium.iteritems() if
                                                                                  k in reconstructed_model.reactions})
                if sol.message == "optimal" and sol.fobj > biomass_target_flux:
                    del gapfilled_reactions[r_id]
                    reconstructed_model.remove_reaction(r_id)
                    removed_reactions.append(r_id)

            if len(removed_reactions) > 0:
                print "Reactions '{}' were removed from list of gapfilled reactions".format(
                    "', '".join(removed_reactions))

            # Renamed gapfilled reactions
            print "Gapfilled {}/{} reactions".format(len(gapfilled_reactions), gapfilled_reactions_len)

            # Update the model for next medium
            models_gf[model_id] = reconstructed_model

            # Write results into tab-delimited file
            with open(output_path, "a") as output_f:
                for r_id, r in sorted(gapfilled_reactions.iteritems(), key=lambda x: x[0]):
                    reaction_annotation = agora_reactions[agora_reactions.ReactionID == r_id]
                    reaction_annotation = "\t".join(reaction_annotation[["ReactionName", "EC", "Subsystem"]].values[0])
                    line = "\t".join([model_id, medium_id, str(medium_no), r.id, r.to_equation_string()])
                    print line
                    output_f.write("{}\t{}\t1\n".format(line, reaction_annotation))
        elif is_growing_insilico and not is_growing_invitro:
            pass


print "Reconstruction completed: {}/{}".format(sum([len(models_f[m].reactions) < len(models_gf[m].reactions) for m in models_f]), len(models_f))