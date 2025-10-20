def get_samples(df, *, included=None, batch=None, condition=None):
    """Get sample names based on inclusion/exclusion, batch and contition.
    
        Args:
            df (pd.DataFrame): DataFrame containing sample information.
            included (bool): If True, return included samples; if False, return excluded samples.
            batch (str or None): If provided, filter samples by this batch.
            condition (str or None): If provided, filter samples by this experimental condition.
        Returns:
            list: List of sample names matching the criteria.
    """
    filtered_df = df.copy()
    
    # Filtro por inclusión
    if included is True:
        filtered_df = filtered_df[filtered_df['included'] == 1]
    elif included is False:
        filtered_df = filtered_df[filtered_df['included'] == 0]
    
    # Filtro por batch
    if batch:
        filtered_df = filtered_df[filtered_df['batch'] == batch]
    
    # Filtro por condition
    if condition:
        filtered_df = filtered_df[filtered_df['condition'] == condition]
    
    samples_subset = filtered_df['sample'].tolist()
    
    print("Samples selected based on criteria:", samples_subset)
    
    return samples_subset

def get_resource(rule,resource) -> int:
	'''
    Retrieve resources for a given rule from config. 	
    '''

	try:
		return config['resources'][rule][resource]
	except KeyError: 
		print(f'Failed to resolve resource for {rule}/{resource}: using default parameters')
		return config["resources"]['default'][resource]
      

# TODO: CHEQUEAR TODO ESTO.
def get_count_table(project, mageck_test):
    """Get the path to the count table based on the normalization/CNV correction method.
    
        Args:
            project (str): Project name.
            mageck_test (str): Type of MAGECK test or normalization method.
        Returns:
            str: Path to the corresponding count table.
    """

    mageck_counts_norms = ["median", "total", "off-target", "non_essen_genes"]

    if mageck_test == "cnvcorr":
        return f"results/cnv_correction/{project}_cnvcorr.count_normalized.txt"
    elif mageck_test in mageck_counts_norms:
        return f"results/mageck_normalize/{mageck_test}/{project}_{mageck_test}.count_normalized.txt"
    else:
        raise ValueError(f"Unknown mageck test option: {mageck_test}")
    
# TODO: AÑADIR EL RESTO DE FUNCIONES