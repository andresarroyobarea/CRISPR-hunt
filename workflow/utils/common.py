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