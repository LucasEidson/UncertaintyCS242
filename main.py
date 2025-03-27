"""Exact Inference: Calculate the exact probability of the query variables given some prior
A query can be answered in a bayes net by computing sums of products of conditional probabilites from the network
"""

""" Rejection Sampling: """
""" Gibbs sampling: """
"""TODO: figure out how to index P(X | Y,Z,... ) - Check Discord, then move on to exact inference"""


def main():
    input = open("easy1.txt").readlines()
    vars = get_bayesian_net(input)
    for v in vars.keys():
        print(v)
        cpt = vars[v].cpt
        for a in cpt:
            for b in a:
                print(b, end=" ")
            print()
    print("Get P(X2 = 1| X1 = 1): ")
    cpt = vars["X2"].cpt
    print(cpt[1][1])
    cpt = vars["X3"].cpt
    print("Get P(X3 = 1| X1 = 1, X2 = 0)")
    print(cpt[2][1])


def get_bayesian_net(lines):
    index = 0
    num_vars = int(lines[index])
    index += 1
    vars = dict()
    for i in range(num_vars):
        var_info = lines[index].split()
        index += 1
        new_var = variable(var_info[0], var_info[1:])
        vars[new_var.name] = new_var
    num_tables = int(lines[index])
    index += 1
    for _ in range(num_tables):
        index += 1  # Skip the blank line in-between each table
        cpt_vars = lines[index].split()
        index += 1
        num_cpt_vars = len(cpt_vars)
        # Transform CPT Vars to variable obj instead of strings:
        for i in range(num_cpt_vars):
            cpt_vars[i] = vars[cpt_vars[i]]
        # first var is child of other vars:
        if num_cpt_vars > 1:
            cpt_vars[0].addParents(cpt_vars[1:])
        cpt = list()
        for _ in range(2 ** (num_cpt_vars - 1)):
            vals = lines[index].split()
            index += 1
            for i in range(len(vals)):
                vals[i] = float(vals[i])
            cpt.append(vals)
        cpt_vars[0].cpt = cpt
    return vars


class variable:
    """Represents a random variable in the bayesian net"""

    def __init__(self, name, domain):
        self.name = name
        self.domain = domain
        self.parents = list()
        self.children = list()
        self.cpt = list()

    def addParents(self, new_parents):
        """Takes a list of variables and adds them to this variables parent list"""
        self.parents = self.parents.extend(new_parents)


if __name__ == "__main__":
    main()
