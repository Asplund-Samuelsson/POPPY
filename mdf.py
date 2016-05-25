#!/usr/bin/env python3

# gibbr commit: bd01559

# Import modules
from scipy import optimize
import numpy as np
import sys, os
import pandas as pd
import collections
import itertools
import argparse

# Import scripts
sys.path.append("/home/johannes/proj/main/tools/pykegg/")
import pykegg

# Define functions
def sWrite(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def sError(string):
    sys.stderr.write(string)
    sys.stderr.flush()


def read_reactions(reactions_text):
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
    # Return stoichiometric matrix (pd.DataFrame)
    return pd.DataFrame(d).fillna(0).astype(int)

def test_read_reactions():
    from pandas.util.testing import assert_frame_equal
    reactions_text = """R1\tA + B <=> C + 2 D\nR2\tA + E <=> 2 F\nR3\tE <=> B + G\n"""
    compounds = ['A','B','C','D','E','F','G']
    exp_df = pd.DataFrame(
        [
            [-1, -1,  0],
            [-1,  0,  1],
            [ 1,  0,  0],
            [ 2,  0,  0],
            [ 0, -1, -1],
            [ 0,  2,  0],
            [ 0,  0,  1]
        ],
        columns = ['R1','R2','R3'],
        index = compounds
    )
    pd.util.testing.assert_frame_equal(exp_df, read_reactions(reactions_text))


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

def test_read_reaction_drGs():
    t1 = "\n".join([
        "R1\t4.0\t-10.0",
        "R1\t5.0\t-9.0",
        "R1\t6.0\t-7.8",
        "R2\t4.0\t2.1",
        "R2\t5.0\t2.5",
        "R2\t6.0\t3.2",
        "R3\t4.0\t10.5",
        "R3\t5.0\t14.6",
        "R3\t6.0\t18.2"
        ]) + "\n"
    exp1 = pd.DataFrame(
        [
            ["R1","4.0",-10.0],
            ["R1","5.0",-9.0],
            ["R1","6.0",-7.8],
            ["R2","4.0",2.1],
            ["R2","5.0",2.5],
            ["R2","6.0",3.2],
            ["R3","4.0",10.5],
            ["R3","5.0",14.6],
            ["R3","6.0",18.2],
        ],
        columns = ["rxn_id","v0","drG"]
    )
    pd.util.testing.assert_frame_equal(exp1, read_reaction_drGs(t1))
    t2 = "R1\ta\t4.0\t23.6\nR2\tb\t5.0\t26.8\n"
    exp2 = pd.DataFrame(
        [
            ["R1",'a',"4.0",23.6],
            ["R2",'b',"5.0",26.8]
        ],
        columns = ["rxn_id","v0","v1","drG"]
    )
    pd.util.testing.assert_frame_equal(exp2, read_reaction_drGs(t2))


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

def test_read_constraints():
    # Test 1: Regular constraints, no range, no step number
    con_text = "C00001\t1e-7\t1e-6\nC00002\t1e-7\t0.001\nC00003\t1e-7\t1\nX00004\t0.5\t0.5\n"
    exp_df = pd.DataFrame(
        [
            ["C00001",0.0000001,0.000001, None, 'lin'],
            ["C00002",0.0000001,0.001, None, 'lin'],
            ["C00003",0.0000001,1.0, None, 'lin'],
            ["X00004",0.5,0.5, None, 'lin']
        ],
        columns = ["cpd_id","x_min","x_max", "steps", "spacing"]
    )
    pd.util.testing.assert_frame_equal(exp_df, read_constraints(con_text))
    # Test 2: Mix of constraints - regular and those with steps
    con_text = "\n".join([
        "A\t0.0001\t0.01",
        "B\t0.0005\t0.5\t5\tlog",
        "C\t0.01\t1\t100\tlin"
    ])
    exp_df = pd.DataFrame(
        [
            ["A",0.0001,0.01, None, 'lin'],
            ["B",0.0005,0.5, 5, 'log'],
            ["C",0.01,1.0, 100, 'lin'],
        ],
        columns = ["cpd_id","x_min","x_max", "steps", "spacing"]
    )
    pd.util.testing.assert_frame_equal(exp_df, read_constraints(con_text))


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

def test_read_ratio_constraints():
    # Test 1; Fixed ratios
    ratio_text = "\n".join([
        "C1\tC2\t10",
        "C3\tC4\t0.1"
    ])
    exp_df = pd.DataFrame(
        [
            ["C1","C2",10,None,5, 'lin'],
            ["C3","C4",0.1,None,5, 'lin']
        ],
        columns = ["cpd_id_num","cpd_id_den","ratio",
                   "ratio_upper", "ratio_step", "spacing"]
    )
    pd.util.testing.assert_frame_equal(exp_df, read_ratio_constraints(ratio_text))
    # Test 2; Ratio ranges
    ratio_text = "\n".join([
        "C1\tC2\t10\t100\t6\tlin",
        "C3\tC4\t0.01\t1\t40\tlog",
        "C5\tC6\t0.1\t1",
        "C7\tC8\t5"
    ])
    exp_df = pd.DataFrame(
        [
            ["C1","C2",10,100,6,'lin'],
            ["C3","C4",0.01,1,40,'log'],
            ["C5","C6",0.1,1,5,'lin'],
            ["C7","C8",5.0,None,5,'lin']
        ],
        columns = ["cpd_id_num","cpd_id_den","ratio",
                   "ratio_upper", "ratio_step", "spacing"]
    )
    pd.util.testing.assert_frame_equal(exp_df, read_ratio_constraints(ratio_text))


def mdf_c(S):
    """Constructs the MDF c vector."""
    # c is all zeroes except for the final element, which represents the MDF
    return np.array([0]*S.shape[0] + [1])

def test_mdf_c():
    S1 = pd.DataFrame([[1,0,-1],[-1,1,-1],[0,-1,2]]) # m=3 cpds, n=3 rxns
    S2 = pd.DataFrame([[1,0,0,-1],[0,-1,2,0]]) # m=2 cpds, n=4 rxns
    np.testing.assert_array_equal(np.array([0,0,0,1]), mdf_c(S1))
    np.testing.assert_array_equal(np.array([0,0,1]), mdf_c(S2))


def mdf_A(S):
    """Constructs the MDF A matrix."""
    # Transposed S in top left corner
    A = np.array(S.T)
    # Below; pos. and neg. identity matrices with same row number as S
    A = np.concatenate((A, np.eye(S.shape[0]), -np.eye(S.shape[0])), axis=0)
    # Pad right side with 1's (S column number) and 0's (2 x S row number)
    A = np.column_stack((A, np.array([1]*S.shape[1] + [0]*S.shape[0]*2)))
    return np.matrix(A)

def test_mdf_A():
    S1 = pd.DataFrame([[1,0,-1],[-1,1,-1],[0,-1,2]]) # m=3 cpds, n=3 rxns
    S2 = pd.DataFrame([[1,0,0,-1],[0,-1,2,0]]) # m=2 cpds, n=4 rxns
    A1 = np.matrix([
        [ 1, -1,  0, 1],
        [ 0,  1, -1, 1],
        [-1, -1,  2, 1],
        [ 1,  0,  0, 0],
        [ 0,  1,  0, 0],
        [ 0,  0,  1, 0],
        [-1,  0,  0, 0],
        [ 0, -1,  0, 0],
        [ 0,  0, -1, 0]
    ])
    A2 = np.matrix([
        [ 1,  0,  1],
        [ 0, -1,  1],
        [ 0,  2,  1],
        [-1,  0,  1],
        [ 1,  0,  0],
        [ 0,  1,  0],
        [-1,  0,  0],
        [ 0, -1,  0]
    ])
    np.testing.assert_array_equal(A1, mdf_A(S1))
    np.testing.assert_array_equal(A2, mdf_A(S2))


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

def test_mdf_b():
    RT = 8.31e-3 * 298.15
    S = read_reactions("R1\tA + B <=> 2 C\nR2\tC + D <=> A\nR3\tD <=> A + B\n")
    drGs = read_reaction_drGs("R2\t-1\nR1\t5\nR3\t-5\n") # Deliberate disorder
    constraints = read_constraints("B\t1e-7\t0.001\nA\t0.0005\t0.05\n") # Deliberate disorder
    expected_mdf_b = np.array([
        -5.0/RT, 1.0/RT, 5.0/RT,
        np.log(0.05), np.log(0.001), np.log(0.01), np.log(0.01),
        -np.log(0.0005), -np.log(1e-7), -np.log(1e-6), -np.log(1e-6)
    ])
    np.testing.assert_array_equal(expected_mdf_b, mdf_b(S, drGs, constraints))


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

def test_mdf_A_eq():
    ratio_df = pd.DataFrame(
        [
            ["C1","C2",10],
            ["C3","C4",0.1]
        ],
        columns = ["cpd_id_num","cpd_id_den","ratio"]
    )
    S = read_reactions("\n".join([
    "R1\tC1 + A <=> C2 + B",
    "R2\tC3 + C0 <=> C4 + D",
    "R3\tC4 + E <=> 2 F"
    ]))
    exp_A_eq = np.matrix([
        [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0]
    ])
    np.testing.assert_array_equal(exp_A_eq, mdf_A_eq(S, ratio_df))


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

def test_mdf_b_eq():
    ratio_df = pd.DataFrame(
        [
            ["C1","C2",10],
            ["C3","C4",0.1]
        ],
        columns = ["cpd_id_num","cpd_id_den","ratio"]
    )
    exp_b_eq = np.array([np.log(10), np.log(0.1)])
    np.testing.assert_array_equal(exp_b_eq, mdf_b_eq(ratio_df))


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

def test_mdf():
    # Define the test case (Markus Janasch's reproduction of glycolysis MDF)
    reactions_text = "\n".join([
        "HXK\talpha-D-glucose[c] + ATP[c] <=> ADP[c] + alpha-D-glucose_6-phosphate[c]",
        "PGI\talpha-D-glucose_6-phosphate[c] <=> beta-D-fructofuranose_6-phosphate[c]",
        "PFK\tATP[c] + beta-D-fructofuranose_6-phosphate[c] <=> ADP[c] + beta-D-fructofuranose_1,6-bisphosphate[c]",
        "FBA\tbeta-D-fructofuranose_1,6-bisphosphate[c] <=> D-glyceraldehyde_3-phosphate[c] + glycerone_phosphate[c]",
        "TPI\tglycerone_phosphate[c] <=> D-glyceraldehyde_3-phosphate[c]",
        "TDH\t2 D-glyceraldehyde_3-phosphate[c] + 2 NAD[c] + 2 PI[c] <=> 2 3-phospho-D-glyceroyl-phosphate[c] + 2 NADH[c]",
        "PGK\t2 3-phospho-D-glyceroyl-phosphate[c] + 2 ADP[c] <=> 2 3-phospho-D-glycerate[c] + 2 ATP[c]",
        "GPM\t2 3-phospho-D-glycerate[c] <=> 2 2-phospho-D-glycerate[c]",
        "ENO\t2 2-phospho-D-glycerate[c] <=> 2 phosphoenolpyruvate[c] + 2 H2O[c]",
        "PYK\t2 ADP[c] + 2 phosphoenolpyruvate[c] <=> 2 ATP[c] + 2 pyruvate[c]"
    ]) + "\n"
    drGs_text = "\n".join([
        "HXK\t-19.5", "PGI\t2.5", "PFK\t-15", "FBA\t19.8",
        "TPI\t5.5", "TDH\t15.6", "PGK\t-36.8", "GPM\t8.4",
        "ENO\t-8.2", "PYK\t-55.4"
    ]) + "\n"
    con_text = "\n".join([
        "ADP[c]\t5.00E-04\t5.00E-04",
        "ATP[c]\t5.00E-03\t5.00E-03",
        "NAD[c]\t5.00E-03\t5.00E-03",
        "NADH[c]\t5.00E-04\t5.00E-04",
        "PI[c]\t5.00E-03\t5.00E-03",
        "H2O[c]\t1\t1"
    ]) + "\n"
    # Load stoichiometric matrix, drGs and constraints
    S = read_reactions(reactions_text)
    drGs = read_reaction_drGs(drGs_text)
    constraints = read_constraints(con_text)
    # Calculate MDF inputs
    c = mdf_c(S)
    A = mdf_A(S)
    b = mdf_b(S, drGs, constraints)
    # Perform MDF
    mdf_result = mdf(c, A, b)
    # Check the MDF result
    assert mdf_result.success # # Optimization terminated successfully
    assert mdf_result.status == 0 # Optimization terminated successfully
    assert abs(mdf_result.x[-1]*8.31e-3*298.15 - 0.4417) / 0.4417 < 0.05


def ratio_range(row):
    """Create a linear or logarithmic range based on a ratio DataFrame row"""
    # For fixed ratios
    if not np.isfinite(row['ratio_upper']):
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

def test_ratio_range():
    ratio_constraints_text = "\n".join([
        "X\tY\t1\t3\t3", "W\tZ\t0.1\t0.4\t4",
        "A\tB\t0.025\t40\t5\tlog", "C\tD\t0.8\t1.1\t4\tlin",
        "Q\tR\t0.6\t1.4", "M\tN\t2", "K1\tK2\t1\t10\t5\tlog"]) + "\n"
    rats = read_ratio_constraints(ratio_constraints_text)
    exp_1 = np.array([1,2,3])
    exp_2 = np.array([0.1,0.2,0.3,0.4])
    exp_3 = np.array([0.025,0.158113883008,1.0,6.32455532034,40.0])
    exp_4 = np.array([0.8,0.9,1.0,1.1])
    exp_5 = np.array([0.6,0.8,1.0,1.2,1.4])
    exp_6 = np.array([2])
    exp_7 = np.array([1.0, 1.77827941004, 3.16227766017, 5.6234132519, 10.0])
    np.testing.assert_array_almost_equal(ratio_range(rats.iloc[0,:]), exp_1)
    np.testing.assert_array_almost_equal(ratio_range(rats.iloc[1,:]), exp_2)
    np.testing.assert_array_almost_equal(ratio_range(rats.iloc[2,:]), exp_3)
    np.testing.assert_array_almost_equal(ratio_range(rats.iloc[3,:]), exp_4)
    np.testing.assert_array_almost_equal(ratio_range(rats.iloc[4,:]), exp_5)
    np.testing.assert_array_almost_equal(ratio_range(rats.iloc[5,:]), exp_6)
    np.testing.assert_array_almost_equal(ratio_range(rats.iloc[6,:]), exp_7)


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

def test_con_range():
    constraints_text = "\n".join([
        "X\t0.1\t0.3\t3", "W\t0.1\t0.4\t4",
        "A\t0.000001\t0.1\t6\tlog", "C\t0.8\t1.1\t4\tlin",
        "Q\t0.6\t1.4", "M\t2\t2", "K1\t0.01\t0.1\t4\tlog"]) + "\n"
    cons = read_constraints(constraints_text)
    exp_1 = np.array([0.1,0.2,0.3])
    exp_2 = np.array([0.1,0.2,0.3,0.4])
    exp_3 = np.array([0.000001,0.00001,0.0001,0.001,0.01,0.1])
    exp_4 = np.array([0.8,0.9,1.0,1.1])
    exp_5 = None # Step number not specified; regular bounds
    exp_6 = None # Step number not specified; regular bounds
    exp_7 = np.array([0.01, 0.021544347, 0.046415888, 0.1])
    np.testing.assert_array_almost_equal(con_range(cons.iloc[0,:]), exp_1)
    np.testing.assert_array_almost_equal(con_range(cons.iloc[1,:]), exp_2)
    np.testing.assert_array_almost_equal(con_range(cons.iloc[2,:]), exp_3)
    np.testing.assert_array_almost_equal(con_range(cons.iloc[3,:]), exp_4)
    assert con_range(cons.iloc[4,:]) == exp_5
    assert con_range(cons.iloc[5,:]) == exp_6
    np.testing.assert_array_almost_equal(con_range(cons.iloc[6,:]), exp_7)


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

def test_ratio_iter():
    # The ratio step determines the number of steps
    ratio_constraints_text = "X\tY\t1\t3\t3\nW\tZ\t0.1\t0.4\t4\n"
    ratio_constraints = read_ratio_constraints(ratio_constraints_text)
    labels = ['cpd_id_num','cpd_id_den','ratio']
    expected_rats = [
        pd.DataFrame([['X','Y',x[0]],['W','Z',x[1]]], columns = labels) \
        for x in itertools.product(np.linspace(1,3,3), np.arange(0.1,0.5,0.1))
    ]
    for rats in enumerate(list(ratio_iter(ratio_constraints))):
        pd.util.testing.assert_frame_equal(rats[1], expected_rats[rats[0]])
    # When the range passes 1 a flip should occur
    # This should make the step size is "even" when comparing ratios <1 and >1
    ratio_constraints_text = "X\tY\t0.025\t40\t17\tlog\nW\tZ\t0.8\t1.1\t4\tlin\n"
    ratio_constraints = read_ratio_constraints(ratio_constraints_text)
    expected_rats = [
        pd.DataFrame([['X','Y',x[0]],['W','Z',x[1]]], columns = labels) \
        for x in itertools.product(
            [0.025,
             0.0396458293784,
             0.0628716714841,
             0.0997039824159,
             0.158113883008,
             0.250742241125,
             0.397635364384,
             0.630583352447,
             1,
             1.58583317514,
             2.51486685937,
             3.98815929664,
             6.32455532034,
             10.029689645,
             15.9054145753,
             25.2233340979,
             40.0 ],
            [0.8, 0.9, 1.0, 1.1]
            )
    ]
    for rats in enumerate(list(ratio_iter(ratio_constraints))):
        pd.util.testing.assert_frame_equal(rats[1], expected_rats[rats[0]])


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

def test_con_iter():
    con_text = "\n".join([
        "A\t0.1\t0.3\t3\tlin",
        "B\t0.01\t1\t3\tlog",
        "C\t0.1\t0.5"
    ])
    constraints = read_constraints(con_text)
    A_lines = ["A\t0.1\t0.1", "A\t0.2\t0.2", "A\t0.3\t0.3"]
    B_lines = ["B\t0.01\t0.01", "B\t0.1\t0.1", "B\t1.0\t1.0"]
    C_lines = ["C\t0.1\t0.5"]
    expected_cons = [
        read_constraints("\n".join(lines)).iloc[:,:3] for lines in \
        itertools.product(A_lines, B_lines, C_lines)
    ]
    multiple_constraints = list(con_iter(constraints))
    for n in range(len(expected_cons)):
        pd.util.testing.assert_frame_equal(multiple_constraints[n], expected_cons[n])


def calc_drGs(S, drGs_std, log_conc, T=298.15, R=8.31e-3):
    """Calculate reaction Gibbs energies"""
    drGs = []
    for rxn_id in list(S.columns):
        drG = float(drGs_std[drGs_std.rxn_id == rxn_id]['drG']) \
        + float(sum(S[rxn_id].T * log_conc)*T*R)
        drGs.append(drG)
    return drGs

def test_calc_drGs():
    S = read_reactions("R1\tA + B <=> C\nR2\tC <=> 2 D\n")
    drGs_std = read_reaction_drGs("R1\t5\nR2\t-2\n")
    log_conc = [-1,-2,-3,-4]
    assert calc_drGs(S, drGs_std, log_conc) == [
        5-(-1-2--3)*298.15*8.31e-3,
        -2-(-3-2*-4)*298.15*8.31e-3
    ]


def multi_mdf(S, all_drGs, constraints, ratio_constraints=None,
              all_directions=False, T=298.15, R=8.31e-3):
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
    all_directions : bool, optional
        Set to True to calculate MDF for all possible reaction direction
        combinations. Not recommended for sets of reactions >20.
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
        A = mdf_A(S_mod)
        b = mdf_b(S_mod, drGs, constraints_mod)

        # Use equality (ratio) constraints if they were specified
        if rats is not None:
            A_eq = mdf_A_eq(S_mod, rats)
            b_eq = mdf_b_eq(rats)
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

def test_multi_mdf_1():
    # TEST CASE 1: COMPLEX
    # Test text input
    reactions_text = "\n".join([
        "R1\tA + X <=> B + Y",
        "R2\tB + C <=> D + Z",
        "R3\tD + X <=> E + C + Y"
    ]) + "\n"
    reaction_drGs_text = "\n".join([
        "R1\t7.0\t-15",
        "R1\t8.0\t-25",
        "R2\t7.0\t20",
        "R2\t8.0\t30",
        "R3\t7.0\t-10",
        "R3\t8.0\t-5"
    ])
    constraints_text = "C\t0.0001\t0.002\n"
    ratio_constraints_text = "X\tY\t1\t3\t3\n"

    # Test input variables setup
    S = read_reactions(reactions_text)
    drGs = read_reaction_drGs(reaction_drGs_text)
    cons = read_constraints(constraints_text)
    rats = read_ratio_constraints(ratio_constraints_text)

    # Perform multiple MDF analyses
    multi_mdf_result = multi_mdf(
        S = S,
        all_drGs = drGs,
        constraints = cons,
        ratio_constraints = rats,
        all_directions = True
    )

    # Expected multiple MDF result
    directions = [
        [1,1,1],                         # All forward
        [1,1,-1], [1,-1,1], [-1,1,1],    # One reverse
        [1,-1,-1], [-1,1,-1], [-1,-1,1], # Two reverse
        [-1,-1,-1]                       # All reverse
    ]
    pHs = ['7.0', '8.0']
    ratios_X_Y = [1, 2, 3]
    column_labels = [
        'v0', 'drGstd_R1', 'drGstd_R2', 'drGstd_R3', 'ratio_X_Y',
        'dir_R1', 'dir_R2', 'dir_R3', 'c_A', 'c_B', 'c_C', 'c_D', 'c_E',
        'c_X', 'c_Y', 'c_Z', 'drGopt_R1', 'drGopt_R2', 'drGopt_R3',
        'success', 'MDF'
    ]
    exp_multi_mdf_result = pd.DataFrame(columns=column_labels)
    for params in itertools.product(directions, pHs, ratios_X_Y):
        cons_mod = cons
        # Change directions
        S_mod = S * params[0]
        drGs_mod = drGs.loc[drGs.v0 == params[1]]
        drGs_mod.is_copy = False
        drGs_mod.loc[:,['drG']] = drGs_mod['drG'] * params[0]
        # Create ratio DataFrame
        rats_mod = pd.DataFrame([['X','Y',params[2]]], columns = list(rats.columns[:3]))
        # Set up MDF inputs
        c = mdf_c(S_mod)
        A = mdf_A(S_mod)
        b = mdf_b(S_mod, drGs_mod, cons_mod)
        A_eq = mdf_A_eq(S_mod, rats_mod)
        b_eq = mdf_b_eq(rats_mod)
        # Calculate MDF
        mdf_result = mdf(c, A, b, A_eq, b_eq)
        # Format and append
        mdf_row = [
            params[1], float(drGs_mod[drGs_mod.rxn_id == 'R1']['drG']),
            float(drGs_mod[drGs_mod.rxn_id == 'R2']['drG']),
            float(drGs_mod[drGs_mod.rxn_id == 'R3']['drG']),
            params[2], params[0][0], params[0][1], params[0][2]
        ]
        if mdf_result.success:
            mdf_row.extend([
                *np.exp(mdf_result.x[:-1]), # Concentrations
                *[
                float(sum(
                    S_mod[x].T * mdf_result.x[:-1])*298.15*8.31e-3 \
                    + drGs_mod[drGs_mod.rxn_id == x]['drG']
                    )
                    for x in ['R1','R2','R3']
                ],
                1,
                mdf_result.x[-1]*298.15*8.31e-3
            ])
        else:
            mdf_row.extend([
                *[None]*8,
                *[None]*3,
                0,
                None
            ])
        # Append row to expected result
        exp_multi_mdf_result = exp_multi_mdf_result.append(pd.DataFrame([mdf_row], columns = column_labels))

    sort_labels = [column_labels[x] for x in [0,4,5,6,7]]
    exp_multi_mdf_result.sort_values(sort_labels, inplace=True)

    pd.util.testing.assert_frame_equal(exp_multi_mdf_result, multi_mdf_result)

def test_multi_mdf_2():
    # TEST CASE 2: SIMPLE
    S = read_reactions("R1\tA + B <=> C\nR2\tC <=> 2 D\n")
    drGs = read_reaction_drGs("R1\t5\nR2\t-2\n")
    constraints = read_constraints("")
    c = mdf_c(S)
    A = mdf_A(S)
    b = mdf_b(S, drGs, constraints)
    mdf_result = mdf(c, A, b)
    exp_mdf_row = [
        5.0, -2.0, 1, 1, *np.exp(mdf_result.x[:-1]),
        float(sum(S['R1'].T * mdf_result.x[:-1])*298.15*8.31e-3 + 5),
        float(sum(S['R2'].T * mdf_result.x[:-1])*298.15*8.31e-3 - 2),
        1, mdf_result.x[-1]*298.15*8.31e-3
    ]
    multi_mdf_result = multi_mdf(S, drGs, constraints)
    assert multi_mdf_result.shape == (1, 12)
    acquired_mdf_row = list(multi_mdf_result.iloc[0,:])
    np.testing.assert_almost_equal(acquired_mdf_row, exp_mdf_row)

def test_multi_mdf_3():
    # TEST CASE 3: RANGE OF FIXED CONCENTRATIONS
    S = read_reactions("R1\tA + B <=> C\nR2\tC <=> 2 D\n")
    drGs = read_reaction_drGs("R1\tX\tZ\t5\nR2\tX\tZ\t-2\n")
    constraints = read_constraints("A\t0.001\t0.003\t3\nC\t0.001\t0.01\t2\tlog")
    c = mdf_c(S)
    A = mdf_A(S)
    exp_mdf_rows = []
    for A_C_con in itertools.product([0.001,0.002,0.003],[0.001,0.01]):
        cA = A_C_con[0]
        cC = A_C_con[1]
        cons = read_constraints("A\t%s\t%s\nC\t%s\t%s\n" % (cA, cA, cC, cC))
        b = mdf_b(S, drGs, cons)
        mdf_result = mdf(c, A, b)
        exp_mdf_row = [
            'X', 'Z', 5.0, -2.0, 1.0, 1.0,
            cA, np.exp(mdf_result.x[1]), cC, np.exp(mdf_result.x[3]),
            float(sum(S['R1'].T * mdf_result.x[:-1])*298.15*8.31e-3 + 5),
            float(sum(S['R2'].T * mdf_result.x[:-1])*298.15*8.31e-3 - 2),
            1.0,
            mdf_result.x[-1]*298.15*8.31e-3
        ]
        exp_mdf_rows.append(exp_mdf_row)
    exp_mdf_result = pd.DataFrame(exp_mdf_rows,
        columns = [
            'v0','v1',
            'drGstd_R1','drGstd_R2',
            'dir_R1','dir_R2',
            'c_A','c_B','c_C','c_D',
            'drGopt_R1','drGopt_R2',
            'success','MDF'
        ],
        index = [0] * len(exp_mdf_rows)
    )
    multi_mdf_result = multi_mdf(S, drGs, constraints)
    assert multi_mdf_result.shape == (6, 14)
    print("")
    print(exp_mdf_result)
    print(multi_mdf_result)
    pd.util.testing.assert_frame_equal(exp_mdf_result, multi_mdf_result)


# Main code block
def main(reaction_file, std_drG_file, outfile_name, cons_file, ratio_cons_file,
         all_directions, T=298.15, R=8.31e-3):

    # Load stoichiometric matrix
    sWrite("\nLoading stoichiometric matrix...")
    S = read_reactions(open(reaction_file, 'r').read())
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

    sWrite("Performing MDF optimization...")
    mdf_table = multi_mdf(S, std_drGs, constraints, ratio_constraints,
                          all_directions, T, R)
    sWrite("\n")

    # Write MDF results to outfile
    sWrite("Saving MDF results to csv...")
    mdf_table.to_csv(outfile_name, na_rep='NA', index=False, float_format='%.10f')
    sWrite(" Done.\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('reactions', type=str,
                        help='Load reactions.')
    parser.add_argument('std_drG', type=str,
                        help='Load standard reaction Gibbs energies.')
    parser.add_argument('outfile', type=str,
                        help='Write MDF table in csv format.')
    parser.add_argument('--constraints', type=str,
                        help='Load concentration bound constraints.')
    parser.add_argument('--ratios', type=str,
                        help='Load concentration ratio constraints.')
    parser.add_argument('--all_directions', action='store_true',
                        help='Analyze MDF for all reaction directions.')
    parser.add_argument('-T', type=float, default=298.15,
                        help='Temperature (K).')
    parser.add_argument('-R', type=float, default=8.31e-3,
                        help='Universal gas constant (kJ/(mol*K)).')
    args = parser.parse_args()
    main(
        args.reactions, args.std_drG, args.outfile, args.constraints,
        args.ratios, args.all_directions, args.T, args.R
    )
