#!/usr/bin/env python3

# Import modules
from scipy import optimize
import numpy as np
import sys, os
import pandas as pd
import collections
import itertools
import argparse

# Import scripts
import pykegg

# Define functions
def sWrite(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def sError(string):
    sys.stderr.write(string)
    sys.stderr.flush()


def read_reactions(reactions_text, proton_name = "C00080"):
    """Create stoichiometric matrix DataFrame from text

    ARGUMENTS

    reactions_text : string
        Tab-delimited string with reaction IDs in the first column
        and KEGG style reaction strings in the second column.

    RETURNS

    pandas.DataFrame
        Pandas DataFrame that corresponds to the stoichiometric matrix. Column
        names are reaction IDs and row indices are compound names.
    """
    # Iterate over lines in reactions_text
    d = collections.OrderedDict()
    for line in reactions_text.split("\n"):
        if line == "":
            continue
        # Extract reaction ID and parsed equation from line
        line = line.rstrip().split("\t")
        rxn_id = line[0]
        eq = pykegg.parse_equation(line[1])
        # List the coefficients
        coefficients = [c[0]*-1 for c in eq[0]] + [c[0] for c in eq[1]]
        # List the compound IDs
        compounds = [c[1] for c in eq[0]] + [c[1] for c in eq[1]]
        # Store a Pandas series for the reaction
        d[rxn_id] = pd.Series(coefficients, index = compounds)
    # Construct stoichiometric matrix (pd.DataFrame)
    S = pd.DataFrame(d).fillna(0).astype(int)
    # Remove protons
    return S.loc[S.index != proton_name]


def read_reaction_drGs(reaction_drGs_text):
    """Create reaction standard Gibbs energies DataFrame from text

    ARGUMENTS

    reaction_drGs_text : string
        Tab-delimited string with reaction IDs in the first column, optional
        values for grouping in intermediate columns and reaction standard Gibbs
        energies in the last column.

    RETURNS

    pandas.DataFrame
        Pandas DataFrame with one column per column in the input argument. The
        first columns contain strings whereas the last column contains the
        reaction standard Gibbs energies in float format.
    """
    # Iterate over lines in reaction_drGs_text
    d = []
    cols = []
    for line in reaction_drGs_text.split("\n"):
        if line == "":
            continue
        # Extract reaction ID, drG and intermediate identifier values
        line = line.rstrip().split("\t")
        rxn_id = line[0]
        drG = line[-1]
        vals = []
        # Set up the column name list
        if cols == []:
            cols = ["rxn_id"]\
            + ["v" + str(i) for i in range(len(line[1:-1]))]\
            + ["drG"]
        for val in line[1:-1]:
            vals.append(val)
        # Store a Pandas series for the line
        d.append([rxn_id] + vals + [drG])
    # Return drG pd.DataFrame
    df = pd.DataFrame(d, columns=cols)
    df["drG"] = pd.to_numeric(df["drG"])
    return df


def read_constraints(constraints_text, default_spacing='lin'):
    """Create constraints DataFrame from text

    ARGUMENTS

    constraints_text : string
        Tab-delimited string with compound IDs in the first column, lower
        concentration bounds in the second column, in M, and upper concentration
        bounds in the third column, in M. Optional fourth and fifth columns
        specify the number of steps and spacing type ('lin' for linear, 'log'
        for logarithmic), respectively.

    RETURNS

    pandas.DataFrame
        Pandas DataFrame with a compound ID column (string), a lower
        concentration bound column (float) and an upper concentration bound
        colunn (float). The fourth and fifth columns contain
    """
    # Iterate over lines in constraints_text
    data = []
    for line in constraints_text.split("\n"):
        row = []
        if line == "":
            continue
        # Extract compound ID and values from line
        line = line.rstrip().split("\t")
        row.append(line[0]) # cpd_id
        row.append(float(line[1])) # x_min
        row.append(float(line[2])) # x_max
        try:
            row.append(int(line[3])) # Ratio range step number
        except IndexError:
            row.append(None) # If step size is missing, use None
        try:
            row.append(str(line[4])) # Concentration range spacing (lin or log)
        except IndexError:
            row.append(str(default_spacing))
        data.append(row)
    # Return constraints (pd.DataFrame)
    return pd.DataFrame(data, columns = [
        "cpd_id", "x_min", "x_max", "steps", "spacing"
    ])


def read_ratio_constraints(ratio_text, default_step=5, default_spacing='lin'):
    """Create ratio constraints DataFrame from text

    ARGUMENTS

    ratio_txt : string
        Tab-delimited string with compound IDs in the first and second column,
        and the ratio of their concentrations (M) in the third column.
        Optionally, the consecutive columns may indicate the upper limit of
        ratios (with column 3 indicating the lower limit), the number of ratios
        and the type of ratio spacing ('lin' for linear spacing and 'log' for
        log10 spacing).

    default_step : int, optional
        If only the lower and upper limits are given for a range of ratios, use
        this value for the step number.

    default_spacing : string, optional
        If no spacing type is specified in ratio_txt, use this spacing value.

    RETURNS

    pandas.DataFrame
        Pandas DataFrame with two compound ID columns (string), a lower limit
        concentration ratio column (float), an upper limit concentration ratio
        column (float) and the concentration ratio range step number (int). The
        third column is interpreted as the fixed ratio when the fourth column
        contains a None value. The last column indicates the type of spacing to
        use for ratio ranges (linear or logarithmic).
    """
    # Iterate over lines in constraints_text
    d = []
    for line in ratio_text.split("\n"):
        if line == "":
            continue
        # Extract compound ID and values from line
        R = []
        line = line.rstrip().split("\t")
        R.append(line[0]) # Numerator
        R.append(line[1]) # Denominator
        R.append(float(line[2])) # Ratio (or lower ratio limit)
        try:
            R.append(float(line[3])) # Upper ratio limit
        except IndexError:
            R.append(None) # If step size is missing, use None
        try:
            R.append(int(line[4])) # Ratio range step number
        except IndexError:
            R.append(int(default_step)) # If step size is missing, use default
        try:
            R.append(str(line[5])) # Ratio range spacing (lin or log)
        except IndexError:
            R.append(str(default_spacing))
        d.append(R)
    # Return stoichiometric matrix (pd.DataFrame)
    return pd.DataFrame(d, columns = ["cpd_id_num","cpd_id_den","ratio",
                        "ratio_upper", "ratio_step", "spacing"])


def mdf_c(S):
    """Constructs the MDF c vector."""
    # c is all zeroes except for the final element, which represents the MDF
    return np.array([0]*S.shape[0] + [1])


def mdf_A(S, net_rxns = []):
    """Constructs the MDF A matrix."""
    # Transposed S in top left corner
    A = np.array(S.T)
    # Below; pos. and neg. identity matrices with same row number as S
    A = np.concatenate((A, np.eye(S.shape[0]), -np.eye(S.shape[0])), axis=0)
    # Pad right side differently depending on whether a network is specified
    if not net_rxns:
        # Pad with 1's (S column number) and 0's (2 x S row number)
        A = np.column_stack((A, np.array([1]*S.shape[1] + [0]*S.shape[0]*2)))
    else:
        # Pad with 1's or 0's (depending on reaction status; S column number),
        # and 0's (2 x S row number)
        mdf_vector = [0 if R in net_rxns else 1 for R in S.columns]
        A = np.column_stack((A, np.array(mdf_vector + [0]*S.shape[0]*2)))
    return np.matrix(A)


def mdf_b(S, drGs, constraints, x_max_default=0.01, x_min_default=0.000001, T=298.15, R=8.31e-3):
    """Constructs the MDF b vector."""
    # Use the stoichiometric matrix to find the order of compounds and reactions
    cpd_order = S.index
    rxn_order = S.columns
    # Initialize numpy array (vector)
    d = np.array([])
    # Add -drG values
    for rxn in rxn_order:
        d = np.append(d, -drGs[drGs.rxn_id == rxn].drG.get_values()[0])
    # Transform drG values to units of RT
    d = d / (R*T)
    # Add x_max values
    for cpd in cpd_order:
        x_max = constraints[constraints.cpd_id == cpd].x_max.get_values()
        if len(x_max) == 0:
            # Use default x_max
            d = np.append(d, np.log(x_max_default))
        else:
            # Use provided x_max
            d = np.append(d, np.log(x_max[0]))
    # Add x_min values
    for cpd in cpd_order:
        x_min = constraints[constraints.cpd_id == cpd].x_min.get_values()
        if len(x_min) == 0:
            # Use default x_max
            d = np.append(d, -np.log(x_min_default))
        else:
            # Use provided x_max
            d = np.append(d, -np.log(x_min[0]))
    return d


def mdf_A_eq(S, ratio_constraints):
    """Construct equality constraints matrix

    ARGUMENTS

    S : pandas.DataFrame
        Pandas DataFrame that corresponds to the stoichiometric matrix. Column
        names are reaction IDs and row indices are compound names.

    ratio_constraints : pandas.DataFrame
        Pandas DataFrame with two compound ID columns (string) and a
        concentration ratio column (float).

    RETURNS

    numpy.matrix
        Equality constraints matrix for concentration ratios. Gives the natural
        logarithm of the ratio between two compounds when multiplied by the
        vector of concentrations (x).
    """
    # Iterate over ratio constraints
    d = []
    for row in ratio_constraints.index:
        a = np.zeros(S.shape[0]+1) # Plus one because of the B (MDF) variable
        np.put(a, S.index.get_loc(ratio_constraints.loc[row, "cpd_id_num"]), 1)
        np.put(a, S.index.get_loc(ratio_constraints.loc[row, "cpd_id_den"]), -1)
        d.append(a)
    return np.matrix(d)


def mdf_b_eq(ratio_constraints):
    """Construct equality constraints vector

    ARGUMENTS

    ratio_constraints : pandas.DataFrame
        Pandas DataFrame with two compound ID columns (string) and a
        concentration ratio column (float).

    RETURNS

    numpy.array
        Equality constraints vector corresponding to the natural logarithms of
        the ratios between compounds specified by ratio_constraints.
    """
    # Iterate over ratio constraints
    b_eq = np.array([])
    for row in ratio_constraints.index:
        b_eq = np.append(b_eq, np.log(ratio_constraints.loc[row, "ratio"]))
    return b_eq


def mdf(c, A, b, A_eq = None, b_eq = None):
    """Perform MDF optimization using the simplex algorithm

    ARGUMENTS

    c : numpy.array
        The linear programming standard form c vector, which for MDF consists of
        zeroes and a single value of one in the last position.
    A : numpy.matrix
        The linear programming standard form A matrix, which for MDF consists of
        a matrix expanded from S in function mdf_A.
    b : numpy.array
        The linear programming standard form b vector, which for MDF consists of
        the standard condition reaction driving forces, the natural logarithms
        of the upper concentration bounds and the negative of the natural
        logarithms of the lower concentration bounds.
    A_eq : numpy.matrix, optional
        Equality constraints matrix for concentration ratios. Gives the natural
        logarithm of the ratio between two compounds when multiplied by the
        vector of concentrations (x).
    b_eq : numpy.array, optional
        Equality constraints vector corresponding to the natural logarithms of
        the ratios between compounds specified by A_eq.

    RETURNS

    scipy.optimize.OptimizeResult
        x : numpy.ndarray
            Vector of optimized concentration natural logarithms. The last value
            of the vector is the MDF in units of RT. Multiply with RT to get the
            value in kJ/mol.
    """
    return optimize.linprog(-c, A_ub=A, b_ub=b, A_eq=A_eq, b_eq=b_eq,
                                  bounds=(None,None))


def ratio_range(row):
    """Create a linear or logarithmic range based on a ratio DataFrame row"""
    # For fixed ratios
    fixed_ratio = False

    try:
        if not np.isfinite(row['ratio_upper']):
            fixed_ratio = True
    except TypeError:
        if row['ratio_upper'] is None:
            fixed_ratio = True

    if fixed_ratio:
        return np.array([row['ratio']])

    # For linear ratio ranges
    if row['spacing'] == 'lin':
        return np.linspace(row['ratio'], row['ratio_upper'], row['ratio_step'])

    # For logarithmic ratio ranges
    if row['spacing'] == 'log':

        # If the range starts or ends at 1.0
        if 1 <= row['ratio'] or 1 >= row['ratio_upper']:
            return np.logspace(np.log10(row['ratio']),
                               np.log10(row['ratio_upper']),
                               row['ratio_step'])

        # If the range is not evenly distributed around one,
        # or has an even step number
        elif 1/row['ratio'] != row['ratio_upper'] or row['ratio_step'] % 2 == 0:
            return np.logspace(np.log10(row['ratio']),
                               np.log10(row['ratio_upper']),
                               row['ratio_step'])

        # The final option is an even distribution and an odd number of steps
        else:
            n_b = int(row['ratio_step'] / 2)
            n_u = row['ratio_step'] - n_b
            b = np.logspace(np.log10(row['ratio']), 0.0, n_b, endpoint = False)
            u = np.logspace(0.0, np.log10(row['ratio_upper']), n_u)
            return np.append(b, u)


def con_range(row):
    """Create a linear or logarithmic range based on a constraints row"""
    # For fixed concentration bounds
    steps = row['steps']
    if steps is None:
        steps = np.nan
    if not np.isfinite(steps):
        return None
    # For linear concentration ranges
    if row['spacing'] == 'lin':
        return np.linspace(row['x_min'], row['x_max'], row['steps'])
    # For logarithmic concentration ranges
    if row['spacing'] == 'log':
        return np.logspace(np.log10(row['x_min']), np.log10(row['x_max']),
                           row['steps'])


def ratio_iter(ratio_constraints):
    """Iterator for constructing DataFrames expected by mdf_A_eq and mdf_b_eq"""
    # Construct row iterator
    row_iter = ratio_constraints.iterrows()

    # Construct ratios column iterator (all possible ratio combinations)
    ratio_col_iter = itertools.product(
        *[ratio_range(row[1]) for row in row_iter]
    )

    # For each possible ratio combination, yield a ratio_constraints DataFrame
    # that is expected by mdf_A_eq and mdf_b_eq
    for ratio_col in ratio_col_iter:
        rats = pd.DataFrame(
            [
                list(ratio_constraints['cpd_id_num']),
                list(ratio_constraints['cpd_id_den']),
                list(ratio_col)
            ],
            index = list(ratio_constraints.columns[:3])
        ).T
        rats.ratio = rats.ratio.astype('float')
        yield rats


def con_iter(constraints):
    """Iterator for constructing DataFrames expected by mdf_A and mdf_b"""

    # Go through the rows and identify those that should not have a range
    non_range_rows = []
    range_rows = []
    for row_n in range(constraints.shape[0]):
        row = constraints.iloc[row_n,]
        if con_range(row) is None:
            non_range_rows.append(row_n)
        else:
            range_rows.append(row_n)

    # Construct concentration column iterator (all possible conc. combinations)
    col_iter = itertools.product(
        *[con_range(row[1]) for row in constraints.iloc[range_rows,].iterrows()]
    )

    # For each possible conc. combination, yield a constraints DataFrame
    # that is expected by mdf_A and mdf_b
    for col in col_iter:
        cons = pd.DataFrame(
            [
                list(constraints.iloc[range_rows,]['cpd_id']),
                list(col), # x_min == x_max
                list(col)  # x_max == x_min
            ],
            index = list(constraints.columns[:3])
        ).T
        cons = cons.append(constraints.iloc[non_range_rows,:3])
        cons.x_min = cons.x_min.astype('float')
        cons.x_max = cons.x_max.astype('float')
        yield cons


def calc_drGs(S, drGs_std, log_conc, T=298.15, R=8.31e-3):
    """Calculate reaction Gibbs energies"""
    drGs = []
    for rxn_id in list(S.columns):
        drG = float(drGs_std[drGs_std.rxn_id == rxn_id]['drG']) \
        + float(sum(S[rxn_id].T * log_conc)*T*R)
        drGs.append(drG)
    return drGs


def multi_mdf(S, all_drGs, constraints, ratio_constraints=None, net_rxns=[],
              all_directions=False, x_max=0.01, x_min=0.000001,
              T=298.15, R=8.31e-3):
    """Run MDF optimization for all condition combinations

    ARGUMENTS

    S : pandas.DataFrame
        Pandas DataFrame that corresponds to the stoichiometric matrix. Column
        names are reaction IDs and row indices are compound names.
    all_drGs : pandas.DataFrame
        Pandas DataFrame with reaction IDs in the first column, condition
        identifier strings in the intermediate columns, and reaction standard
        Gibbs energies in float format in the last column.
    constraints : pandas.DataFrame
        Pandas DataFrame with a compound ID column (string), a lower
        concentration bound column (float) and an upper concentration bound
        colunn (float).
    ratio_constraints : pandas.DataFrame, optional
        Pandas DataFrame with two compound ID columns (string), a lower limit
        concentration ratio column (float), an upper limit concentration ratio
        column (float) and the concentration ratio range step number (int). The
        third column is interpreted as the fixed ratio when the fourth column
        contains a None value. The last column indicates the type of spacing to
        use for ratio ranges (linear or logarithmic).
    net_rxns : list of strings
        List with strings referring to the background network reactions for
        network-embedded MDF analysis (NE-MDF). The reactions should be in S.
    all_directions : bool, optional
        Set to True to calculate MDF for all possible reaction direction
        combinations. Not recommended for sets of reactions >20.
    x_max : float
        Maximum default metabolite concentration (M).
    x_min : float
        Minimum default metabolite concentration (M).
    T : float
        Temperature (K).
    R : float
        Universal gas constant (kJ/(mol*K)).

    RETURNS

    mdf_table : pandas.DataFrame
        A Pandas DataFrame containing all MDF results for a single pathway. Each
        row corresponds to one individual MDF optimization, with the parameters
        described in the columns:
        v0 ... : string
            Condition identifiers as supplied in all_drGs.
        drG_std(rxn_id) : float
            The standard reaction Gibbs energy for the reaction 'rxn_id'.
        [cpd_id_num]/[cpd_id_den] ... : float
            Ratio of concentration between compounds 'cpd_id_num' and
            'cpd_id_den'.
        dir(rxn_id) ... : int
            The direction used for the reaction 'rxn_id'. The order is the same
            as the columns in S.
        [cpd_id] ... : float
            Optimized concentration for compound 'cpd_id' (M).
        drG_opt(rxn_id) : float
            The optimized reaction Gibbs energy for reaction 'rxn_id' (kJ/mol).
        success : int
            Indicates optimization success (1) or failure (0).
        MDF : float
            The Max-min Driving Force determined through linear optimization
            (kJ/mol).
    """
    # All drGs
    # ->  All ratio combinations
    #     ->  All directions

    # Number of reactions
    n_rxn = S.shape[1]

    # List the condition identifiers
    conditions = list(all_drGs.columns[1:-1])

    # Create column labels for output DataFrame
    if ratio_constraints is not None:
        ratio_labels = [
            'ratio_' + ratio_constraints.iloc[row,:]['cpd_id_num'] + \
            '_' + ratio_constraints.iloc[row,:]['cpd_id_den'] \
            for row in range(ratio_constraints.shape[0])
        ]
    else:
        ratio_labels = []

    column_labels = [
        *conditions,
        *['drGstd_' + rxn_id for rxn_id in list(S.columns)],
        *ratio_labels,
        *['dir_' + rxn_id for rxn_id in list(S.columns)],
        *['c_' + cpd_id for cpd_id in list(S.index)],
        *['drGopt_' + rxn_id for rxn_id in list(S.columns)],
        'success',
        'MDF'
    ]

    # Also create labels for sorting (conditions, ratios and directions)
    sort_labels = [
        *conditions,
        *ratio_labels,
        *['dir_' + rxn_id for rxn_id in list(S.columns)]
    ]

    # Iterator preparation
    def prep_iter():

        # Set up conditions iterator
        if len(conditions):
            cond_iter = all_drGs[conditions].drop_duplicates().iterrows()
        else:
            cond_iter = [None]

        # Set up directions iterator
        if not all_directions:
            dir_iter = [[1]*n_rxn]
        else:
            dir_iter = itertools.product([1,-1], repeat=n_rxn)

        # Set up ratios iterator
        if ratio_constraints is not None:
            rats_iter = ratio_iter(ratio_constraints)
        else:
            rats_iter = [None]

        # Set up fixed concentration range constraints iterator
        cons_iter = con_iter(constraints)

        return itertools.product(cond_iter, dir_iter, rats_iter, cons_iter)

    # Set up output DataFrame
    mdf_table = pd.DataFrame(columns = column_labels)

    # Determine number of rows that will be produced
    M = 0
    for i in prep_iter():
        M += 1

    # Iterate over all combinations of conditions, directions and ratios
    n = 0

    for params in prep_iter():
        n += 1
        progress = float(n / M * 100)
        sWrite("\rPerforming MDF optimization... %0.1f%%" % progress)
        # Extract specific condition, direction and ratio constraints
        if params[0] is not None:
            condition = pd.DataFrame(params[0][1]).T
        else:
            condition = None
        direction = params[1]
        rats = params[2]
        constraints_mod = params[3]

        # Obtain specific standard reaction Gibbs energies with correct sign
        if condition is not None:
            drGs = pd.merge(condition, all_drGs)
        else:
            drGs = all_drGs
        drGs.is_copy = False
        drGs.loc[:,['drG']] = drGs['drG'] * direction

        # Modify direction (sign) of reactions in the stoichiometric matrix
        S_mod = S * direction

        # Set up MDF inputs
        c = mdf_c(S_mod)
        A = mdf_A(S_mod, net_rxns)
        b = mdf_b(S_mod, drGs, constraints_mod, x_max, x_min, T, R)

        # Use equality (ratio) constraints if they were specified
        if rats is not None:
            A_eq = mdf_A_eq(S_mod, rats)
            b_eq = mdf_b_eq(rats)
            # If the ratio constraints have been filtered out, set to None
            if not A_eq.size or not b_eq.size:
                A_eq = None
                b_eq = None
        else:
            A_eq = None
            b_eq = None

        # Perform MDF
        mdf_result = mdf(c, A, b, A_eq, b_eq)

        # Prepare conditions list
        if condition is not None:
            conditions_list = list(condition.iloc[0,:])
        else:
            conditions_list = []

        # Prepare ratios list
        if rats is not None:
            rats_list = list(rats.ratio)
        else:
            rats_list = []

        # Format results row
        mdf_row = [
            *conditions_list,
            *[float(drGs[drGs.rxn_id == rxn_id]['drG']) for rxn_id in S_mod.columns],
            *rats_list,
            *direction,
        ]
        if mdf_result.success:
            mdf_row.extend([
                *np.exp(mdf_result.x[:-1]), # Concentrations
                *calc_drGs(S_mod, drGs, mdf_result.x[:-1]), # Reaction Gibbs energies
                1, # Success
                mdf_result.x[-1]*R*T # MDF value
            ])
        else:
            mdf_row.extend([
                *[np.nan]*S_mod.shape[0], # Concentrations
                *[np.nan]*S_mod.shape[1], # Reaction Gibbs energies
                0, # Failure
                np.nan # No MDF value
            ])
        # Append row to expected result
        mdf_table = mdf_table.append(pd.DataFrame([mdf_row], columns = column_labels))

    return mdf_table.sort_values(sort_labels)


# Main code block
def main(reaction_file, std_drG_file, outfile_name, cons_file, ratio_cons_file,
         pw_rxn_file, all_directions, T=298.15, R=8.31e-3, proton_name='C00080',
         x_max_default=0.01, x_min_default=0.000001):

    # Load stoichiometric matrix
    sWrite("\nLoading stoichiometric matrix...")
    S = read_reactions(open(reaction_file, 'r').read(), proton_name)
    sWrite(" Done.\n")

    # Load standard reaction Gibbs energies
    sWrite("Loading standard reaction Gibbs energies...")
    std_drGs = read_reaction_drGs(open(std_drG_file, 'r').read())
    sWrite(" Done.\n")

    # Load constraints
    if cons_file:
        sWrite("Loading concentration bound constraints...")
        constraints = read_constraints(open(cons_file, 'r').read())
        sWrite(" Done.\n")
    else:
        constraints = read_constraints("")
    if ratio_cons_file:
        sWrite("Loading concentration ratio constraints...")
        ratio_constraints = read_ratio_constraints(
            open(ratio_cons_file, 'r').read()
        )
        # Filter the ratio constraints to those relevant to S
        ratio_constraints = ratio_constraints[
                                ratio_constraints['cpd_id_num'].isin(S.index)
                            ]
        ratio_constraints = ratio_constraints[
                                ratio_constraints['cpd_id_den'].isin(S.index)
                            ]
        sWrite(" Done.\n")
    else:
        ratio_constraints = None
    if pw_rxn_file:
        sWrite("Reading pathway reactions for NE-MDF...")
        pw_rxns = list(filter(None,
                       [x.strip() for x in open(pw_rxn_file, 'r').readlines()]))
        # Create list of network reactions rather than pathway reactions
        # for NE-MDF
        net_rxns = list(set(S.columns) - set(pw_rxns))
        sWrite(" Done.\n")
    else:
        net_rxns = []

    sWrite("Performing MDF optimization...")
    mdf_table = multi_mdf(S, std_drGs, constraints, ratio_constraints, net_rxns,
                          all_directions, x_max_default, x_min_default, T, R)
    sWrite("\n")

    # Write MDF results to outfile
    sWrite("Saving MDF results to csv...")
    mdf_table.to_csv(outfile_name, na_rep='NA', index=False, float_format='%.10f')
    sWrite(" Done.\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'reactions', type=str,
        help='Load reactions.'
        )
    parser.add_argument(
        'std_drG', type=str,
        help='Load standard reaction Gibbs energies.'
        )
    parser.add_argument(
        'outfile', type=str,
        help='Write MDF table in csv format.'
        )
    parser.add_argument(
        '--constraints', type=str,
        help='Load concentration bound constraints.'
        )
    parser.add_argument(
        '--ratios', type=str,
        help='Load concentration ratio constraints.'
        )
    parser.add_argument(
        '--pathway', type=str,
        help='Specify pathway reactions for NE-MDF.'
        )
    parser.add_argument(
        '--all_directions', action='store_true',
        help='Analyze MDF for all reaction directions.'
        )
    parser.add_argument(
        '-T', type=float, default=298.15,
        help='Temperature (K).'
        )
    parser.add_argument(
        '-R', type=float, default=8.31e-3,
        help='Universal gas constant (kJ/(mol*K)).'
        )
    parser.add_argument(
        '-H', '--proton_name', default='C00080',
        help='Name used to identify protons.'
        )
    parser.add_argument(
        '--min_conc', type=float, default=0.000001,
        help='Default minimum concentration (M).'
        )
    parser.add_argument(
        '--max_conc', type=float, default=0.01,
        help='Default maximum concentration (M).'
        )
    args = parser.parse_args()
    main(
        args.reactions, args.std_drG, args.outfile, args.constraints,
        args.ratios, args.pathway, args.all_directions, args.T, args.R,
        args.proton_name, args.max_conc, args.min_conc
    )
