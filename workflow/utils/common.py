def get_resource(rule, resource) -> int:
    """
    Function to parse config.yaml to retrieve computational resources for each rule. Returns an int
    """
    try:
        return config["resources"][rule][resource]
    except KeyError:
        print(f"Failed to get resource for {rule}/{resource}: using default parameters")
        return config["resources"]["default"][resource]
    
def get_samples(df, *, included=True, batch=None, condition=None):
    """
    Function to subset samples from the sample.tsv file based on the included, barch and
    condition variables. It could be updated to work with any additional variable.
    """
    
    criteria = []

    if included:
        criteria.append("included == 1")

    if batch:
        criteria.append(f"batch == '{batch}'")

    if condition:
        criteria.append(f"condition == '{condition}'")

    query_str = " and ".join(criteria)
    print(f"Filtering samples with query: {query_str}")

    sample_subset = df.query(query_str)["sample"].unique().tolist()
    print(f"Selected samples: {sample_subset}")

    return sample_subset