from circuit import Circuit, OptimizedCircuit
from matrix import new_vector, new_matrix


# TODO: Evaluate for submission 2A/2B


def generate_fd_analysis_circuit(n: int, to_csv: bool = False, optimized_circuit_class: bool = False) -> (Circuit, int):
    # noinspection GrazieInspection
    """
    Generates a circuit meant to analyse a finite difference mesh of size n.
    Branches within the finite difference matrix are fitted with 10k ohm resistors.
    :param optimized_circuit_class: Optional parameter to use optimized circuit class (banded cholesky algorithms)
    :param n: Size of finite difference mesh
    :param to_csv: Optional flag for generating a csv file of the finite difference mesh. This is technically required
             by the homework spec, but as the Q2 program could (and most likely will) just call this method, I leave
             this as an optional.
    :return: Circuit for analysis of the FD mesh
    """
    n_nodes = n * n
    n_branches = 2 * n * (n - 1) + 1  # Add one for a virtual 'ohm-meter'
    """
    Owing to the limitations of our numerical circuit solver, our ohm-meter needs to
    both conform to the general branch model and have a resistive element.
    We can either use a voltage/resistor series branch or a current/resistor parallel branch.
    I have elected to use the former, with a voltage/resistor pair of 10V/10kΩ. In practice, this produces good 
    results with minimal precision issues.
    """
    J = new_vector(n_branches)
    R = new_vector(n_branches, init_val=10000)  # Branches have 10kΩ resistors
    E = new_vector(n_branches)
    E[int(n_branches / 2)] = 10  # This is the branch that allows us to solve for resistance via voltage analysis

    A = new_matrix(n_nodes, n_branches)

    """
    The (full) incidence matrix is generated by tracing the predictable path that the terminals (positive/negative 1's) 
    leave through the matrix.
    This program uses a node/branch allocation scheme such that path fits a pattern.
    - Nodes are assigned along diagonals following the first norm of its indices in indexed-node notation (i, j).
    - Node 1 is assigned at top left, node n^2 is assigned at bottom right
    - Branches are similarly assigned along diagonals; branches one and two place their positive terminals at (0, 0) and
      the final branches place their negative terminals at (n-1, n-1).
    
    The pattern visible within the terminal path is described thus:
    - Looking between nodes, nonzero elements of a given polarity monotonically increase in branch index. This is to 
      say that if a node i has a +1 in branch j, then no nodes i+k can have a +1 in any branch j-l. This is also true 
      for -1 elements.
    - There are two trends that split the pattern in half. The below statements apply for +1, and match the 
      implementation.
        - Within the first half (by branch), there are two +1 elements per node. This is to say, we can increment
          through the incidence matrix, placing +1 elements two at a time for the first half of all branches. We
          can easily keep track of which row we place elements on via an incrementing pointer.
        - At the "ohm-meter" sentinel branch, we can increment the aforementioned pointer and place a plus at that 
          column/branch.
        - Within the second half, the nodes stagger between one and two +1's per node:
          - This stagger can be implemented by performing the following loop i times for each i from n-1 down to 1, and 
            incrementing the node ptr an extra time between each i:
              - Place a +1
              - Increment the node ptr
              - Place a +1
    - The above process also applies for -1, with the caveat that the iteration occurs backwards and the node pointer
      decrements instead.
      
    The reduced incidence matrix can be gained from the full one by removing the row corresponding to one of the 
    terminals of the ohm-meter branch. The negative one serves nicely, as it is greater than the positive one and
    we need to know the location of the other to calculate resistance of the whole lrn.
    """
    node_plus_ptr = 0
    branch_plus_ptr = 0
    node_minus_ptr = n_nodes - 1
    branch_minus_ptr = n_branches - 1
    for i in range(1, n):
        for j in range(i):
            A[node_plus_ptr][branch_plus_ptr] = 1
            A[node_minus_ptr][branch_minus_ptr] = -1
            branch_plus_ptr += 1
            branch_minus_ptr -= 1
            A[node_plus_ptr][branch_plus_ptr] = 1
            A[node_minus_ptr][branch_minus_ptr] = -1
            branch_plus_ptr += 1
            branch_minus_ptr -= 1
            node_plus_ptr += 1
            node_minus_ptr -= 1

    A[node_plus_ptr][branch_plus_ptr] = -1
    A[node_minus_ptr][branch_minus_ptr] = 1
    pos_node = node_plus_ptr
    neg_node = node_minus_ptr
    branch_plus_ptr += 1
    branch_minus_ptr -= 1
    for i in range(n - 1, 0, -1):
        for j in range(i):
            A[node_plus_ptr][branch_plus_ptr] = 1
            A[node_minus_ptr][branch_minus_ptr] = -1
            branch_plus_ptr += 1
            branch_minus_ptr -= 1
            node_plus_ptr += 1
            node_minus_ptr -= 1
            A[node_plus_ptr][branch_plus_ptr] = 1
            A[node_minus_ptr][branch_minus_ptr] = -1
            branch_plus_ptr += 1
            branch_minus_ptr -= 1

        node_plus_ptr += 1
        node_minus_ptr -= 1

    del A[neg_node]  # Now it should be full rank

    c = OptimizedCircuit(J, R, E, A) if optimized_circuit_class else Circuit(J, R, E, A)
    if to_csv:
        c.load_to_csv(f"incidence{n}.csv")

    return c, pos_node