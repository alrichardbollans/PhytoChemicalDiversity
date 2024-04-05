import os

import pandas as pd
from pkg_resources import resource_filename

from phytochempy.data_compilation_utilities import get_wikidata, get_knapsack_data, merge_and_tidy_compound_datasets, tidy_final_dataset, \
    add_npclassifier_info
from phytochempy.knapsack_searches import get_knapsack_compounds_in_family, tidy_knapsack_results

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')
_tidied_outputs_folder = resource_filename(__name__, 'tidied_outputs')

_output_path = resource_filename(__name__, 'outputs')

all_taxa_compound_csv = os.path.join(_output_path, 'all_taxa_compound_data.csv')

if not os.path.isdir(_temp_outputs_path):
    os.mkdir(_temp_outputs_path)
if not os.path.isdir(_tidied_outputs_folder):
    os.mkdir(_tidied_outputs_folder)
if not os.path.isdir(_output_path):
    os.mkdir(_output_path)

if __name__ == '__main__':
    # Define context
    comp_id_column = 'InChIKey'  # Where appropriate, which column should be used to determine compound uniqueness. This is not applicable to some properties, e.g. where SMILES must be used to generate data
    families = ['Gelsemiaceae', 'Gentianaceae', 'Apocynaceae', 'Loganiaceae', 'Rubiaceae']
    wiki_data_id_for_order = 'Q21754'

    ## Get compound-taxa pair data
    # get_wikidata(wiki_data_id_for_order, os.path.join(_temp_outputs_path, 'wikidata.csv'), os.path.join(_tidied_outputs_folder, 'wikidata.csv'))
    #
    #
    # def _temp_out_for_fam(faml: str) -> str:
    #     return os.path.join(_temp_outputs_path, faml + '_kn_search.csv')
    #
    #
    # def _temp_out_for_fam_Acc(faml: str) -> str:
    #     return os.path.join(_temp_outputs_path, faml + '_kn_search_accepted_info.csv')
    #
    #
    # for fam in families:
    #     get_knapsack_compounds_in_family(fam, _temp_out_for_fam(fam))
    #     tidy_knapsack_results(_temp_out_for_fam(fam), _temp_out_for_fam_Acc(fam), fam, cirpy_cache_dir=_temp_outputs_path,
    #                           add_smiles_and_inchi=True)
    #
    # all_kn_dfs = pd.DataFrame()
    #
    # for fam in families:
    #     new_df = pd.read_csv(_temp_out_for_fam_Acc(fam), index_col=0)
    #     all_kn_dfs = pd.concat([all_kn_dfs, new_df])
    #
    # all_kn_dfs.to_csv(os.path.join(_tidied_outputs_folder, 'knapsack_data.csv'))

    ## Merge and tidy the data
    tidy_wiki_data = pd.read_csv(os.path.join(_tidied_outputs_folder, 'wikidata.csv'), index_col=0)
    tidy_knapsack_data = pd.read_csv(os.path.join(_tidied_outputs_folder, 'knapsack_data.csv'), index_col=0)
    all_compounds_in_taxa = merge_and_tidy_compound_datasets([tidy_wiki_data, tidy_knapsack_data],
                                                             os.path.join(_tidied_outputs_folder, 'merged_data.csv'))

    ## Add extra information related to the compound properties

    # These steps can be included/removed as needed
    # For the longer processes, to avoid repeats you can simply read the associated temp_output if the step has already been run
    with_npclass_classes = add_npclassifier_info(all_compounds_in_taxa, _temp_outputs_path, os.path.join(_tidied_outputs_folder, 'npclassifier.csv'))

    ### Then tidy and output final dataset
    tidy_final_dataset(with_npclass_classes, _tidied_outputs_folder, all_taxa_compound_csv, comp_id_column)

    summary = pd.read_csv(all_taxa_compound_csv, index_col=0)
    summary.describe(include='all').to_csv(os.path.join(_output_path, 'all_taxa_compound_data_summary.csv'))
