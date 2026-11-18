def allocation(stage, key):
    return config["resources"][stage][key]


def quality_input(wildcards):
    return f"{OUT}/checkm2" if config["checkm2"]["exclude_failed"] else []


