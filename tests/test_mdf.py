#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from mdf import *

# Define tests
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


def test_mdf_c():
    S1 = pd.DataFrame([[1,0,-1],[-1,1,-1],[0,-1,2]]) # m=3 cpds, n=3 rxns
    S2 = pd.DataFrame([[1,0,0,-1],[0,-1,2,0]]) # m=2 cpds, n=4 rxns
    np.testing.assert_array_equal(np.array([0,0,0,1]), mdf_c(S1))
    np.testing.assert_array_equal(np.array([0,0,1]), mdf_c(S2))


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

    # Test construction of a network-embedded A matrix, where the MDF is only
    # optimized based on the pathway in focus

    # Name the reactions in S
    S1.columns = ['R1','R2','R3']
    S2.columns = ['R4','R5','R6','R7']

    # Specify the "network"
    net1a = ['R2']
    net1b = ['R1','R3']
    net1c = False
    net2a = ['R5']
    net2b = ['R7','R6']
    net2c = ['R7','R6','R5','R4']

    # Expected A matrices
    A_net1a = np.matrix([
        [ 1, -1,  0, 1],
        [ 0,  1, -1, 0],
        [-1, -1,  2, 1],
        [ 1,  0,  0, 0],
        [ 0,  1,  0, 0],
        [ 0,  0,  1, 0],
        [-1,  0,  0, 0],
        [ 0, -1,  0, 0],
        [ 0,  0, -1, 0]
    ])
    A_net1b = np.matrix([
        [ 1, -1,  0, 0],
        [ 0,  1, -1, 1],
        [-1, -1,  2, 0],
        [ 1,  0,  0, 0],
        [ 0,  1,  0, 0],
        [ 0,  0,  1, 0],
        [-1,  0,  0, 0],
        [ 0, -1,  0, 0],
        [ 0,  0, -1, 0]
    ])
    A_net1c = np.matrix([
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
    A_net2a = np.matrix([
        [ 1,  0,  1],
        [ 0, -1,  0],
        [ 0,  2,  1],
        [-1,  0,  1],
        [ 1,  0,  0],
        [ 0,  1,  0],
        [-1,  0,  0],
        [ 0, -1,  0]
    ])
    A_net2b = np.matrix([
        [ 1,  0,  1],
        [ 0, -1,  1],
        [ 0,  2,  0],
        [-1,  0,  0],
        [ 1,  0,  0],
        [ 0,  1,  0],
        [-1,  0,  0],
        [ 0, -1,  0]
    ])
    A_net2c = np.matrix([
        [ 1,  0,  0],
        [ 0, -1,  0],
        [ 0,  2,  0],
        [-1,  0,  0],
        [ 1,  0,  0],
        [ 0,  1,  0],
        [-1,  0,  0],
        [ 0, -1,  0]
    ])
    np.testing.assert_array_equal(A_net1a, mdf_A(S1, net1a))
    np.testing.assert_array_equal(A_net1b, mdf_A(S1, net1b))
    np.testing.assert_array_equal(A_net1c, mdf_A(S1, net1c))
    np.testing.assert_array_equal(A_net2a, mdf_A(S2, net2a))
    np.testing.assert_array_equal(A_net2b, mdf_A(S2, net2b))
    np.testing.assert_array_equal(A_net2c, mdf_A(S2, net2c))


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


def test_calc_drGs():
    S = read_reactions("R1\tA + B <=> C\nR2\tC <=> 2 D\n")
    drGs_std = read_reaction_drGs("R1\t5\nR2\t-2\n")
    log_conc = [-1,-2,-3,-4]
    assert calc_drGs(S, drGs_std, log_conc) == [
        5-(-1-2--3)*298.15*8.31e-3,
        -2-(-3-2*-4)*298.15*8.31e-3
    ]


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
