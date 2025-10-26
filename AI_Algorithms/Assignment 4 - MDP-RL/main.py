from value_iteration import value_iteration
from q_learning import q_learning
from grid import Grid
from uncertain_grid import Uncertain_Grid

VI_POLICY_0 = {(0, 0): 'east', (0, 1): 'east', (0, 2): 'south', (1, 0): 'north',
                    (1, 1): 'east', (1, 2): 'south', (2, 0): 'north'}

VI_POLICY_1 = {(0, 0): 'east', (0, 1): 'east', (0, 2): 'south', (0, 3): 'south', (0, 4): 'south',
               (1, 0): 'north', (1, 1): 'east', (1, 2): 'south', (1, 3): 'south', (1, 4): 'south',
               (2, 2): 'east', (2, 3): 'south', (2, 4): 'west', (3, 0): 'east', (3, 1): 'south', (3, 2): 'south',
               (3, 3): 'south', (3, 4): 'west', (4, 0): 'east', (4, 1): 'east', (4, 3): 'west', (4, 4): 'west'}


# Test one: Sanity check for value iteration
def sanity_vi(print_policy=False):
    problem_obj = Grid('routes_1.txt')
    vi_policy = value_iteration(problem_obj)
    print("Value iteration has finished running for problem 1, no compilation errors")
    if print_policy:
        problem_obj.pretty_print_policy(vi_policy)
    if vi_policy == VI_POLICY_0:
        print('Correct policy for sanity problem 1')
    else:
        print('Wrong policy for sanity problem 1')
        print('Please check if all functions are completed (marked with TODO comments)')
        print('Exiting sanity check for value iteration')
        return False
    problem_obj = Grid('routes_2.txt', grid_file=1)
    vi_policy = value_iteration(problem_obj)
    print(vi_policy)
    print("Value iteration has finished running for problem 2, no compilation errors")
    if print_policy:
        problem_obj.pretty_print_policy(vi_policy)
    if vi_policy == VI_POLICY_1:
        print('Correct policy for sanity problem 2')
        return True
    else:
        print('Wrong policy for sanity problem 2')
        return False


# Test one: Sanity check for value iteration
def sanity_ql(print_policy=False):
    problem_obj = Uncertain_Grid('routes_1.txt')
    ql_policy = q_learning(problem_obj)
    print("Q-Learning has finished running for problem 1, no compilation errors")
    if print_policy:
        problem_obj.pretty_print_policy(ql_policy)
    if ql_policy == VI_POLICY_0:
        print('Correct policy for sanity problem 1')
    else:
        print('Wrong policy for sanity problem 1')
        print('Please check if all functions are completed (marked with TODO comments)')
        print('Exiting sanity check for value iteration')
        return False
    problem_obj = Uncertain_Grid('routes_2.txt', grid_file=1)
    ql_policy = q_learning(problem_obj, num_episodes=100000)
    print("Q Learning has finished running for problem 2, no compilation errors")
    if print_policy:
        problem_obj.pretty_print_policy(ql_policy)
    if ql_policy == VI_POLICY_1:
        print('Correct policy for sanity problem 2')
        return True
    else:
        print('Wrong policy for sanity problem 2')
        return False


print("\n------ Testing value iteration ------ \n (You may update print_policy=False in order to not print the policy)\n")
vi_success = sanity_vi(print_policy=True)
print("\n------ Testing Q-learning ------ \n (You may update print_policy=False in order to not print the policy)\n")
ql_success = sanity_ql(print_policy=True)

print("\n\n----Summary----")
print("Value iteration success: ", vi_success)
print("Q-learning success: ", ql_success)

