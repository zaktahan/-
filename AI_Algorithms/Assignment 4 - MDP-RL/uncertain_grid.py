import random
import pandas as pd
import numpy as np


class Uncertain_Grid:
    def __init__(self, grid_path, problem_type='south', grid_file=0):
        self.__grid = self.create_grid(grid_path)
        self.start_state, self.end_state = self.create_states_0(grid_file=grid_file)
        self.current_state = self.start_state
        self.__rewards_grid = self.__create_rewards_grid()

    def create_grid(self, path):
        """
        :param path: path to grid file
        :return: grid as matrix
        """
        with open(path) as infile:
            return np.array(np.matrix([list(line.strip()) for line in infile.readlines()]))

    def create_states_0(self, grid_file=0):
        if grid_file == 0:
            start = (0, 0)
            end = (2, 2)
        if grid_file == 1:
            start = (0, 0)
            end = (4, 2)
        if grid_file == 2:
            start = (0, 0)
            end = (6, 3)
        return start, end

    def __create_rewards_grid(self):
        """
        :return: reward of performing an action in state (i,j) for each i,j in grid
        """
        rewards_grid = np.zeros(self.__grid.shape)
        for i in range(self.__grid.shape[0]):
            for j in range(self.__grid.shape[1]):
                if self.__grid[i][j] == 'p':
                    rewards_grid[i][j] = -5
                elif self.__grid[i][j] == 'q':
                    rewards_grid[i][j] = 1
                if (i, j) == self.end_state:
                    rewards_grid[i][j] = 100
        return rewards_grid

    def __get_action_cost(self):
        """
        :return: action cost (hard-coded -1)
        """
        return -1

    def get_reward(self, state, action):
        """
        :param state: (x,y) tuple
        :param action:
        :return: the total reward of performing any action from the given state
        """
        reward = 0
        reward += self.__get_action_cost()
        reward += self.__rewards_grid[state[0]][state[1]]
        return reward

    def __get_actions(self, x, y):
        """
        :param x: x location
        :param y: y location
        :return: nested dictionary in the form {direction:{new_location(x,y):probability}}
        for example: {'north': {(0,0):0.9, (2,0):0.1}, 'south':{(2,0):1.0}}
        """
        # x, y = self.current_state
        pt = {'south': (x + 1, y), 'east': (x, y + 1), 'north': (x - 1, y), 'west': (x, y - 1)}

        # Get valid moves from current location
        final_neighbors = {}
        for direction, location in pt.items():
            if self.__is_valid(location[0], location[1]):
                final_neighbors[direction] = location
            else:
                final_neighbors[direction] = (x, y)  # Stuck in the same place

        # Calculate probabilities for each direction chosen
        directions_probabilities = {}
        for direction, location in final_neighbors.items():
            if direction == 'south':
                directions_probabilities['south'] = {location: 1}
            else:
                if location == final_neighbors['south']:
                    directions_probabilities[direction] = {location: 1}
                else:
                    directions_probabilities[direction] = {location: 0.9}
                    directions_probabilities[direction][final_neighbors['south']] = 0.1
        return directions_probabilities

    def __get_all_states(self):
        """
        :return: list of all the reachable states
        """
        states = self.__create_value_array()
        del states[self.end_state]
        return states.keys()

    def __is_valid(self, x, y):
        """
        :param x: x location
        :param y: y location
        :return: True - if the location is a valid, False - otherwise
        """
        if self.__grid.shape[0] > x >= 0 and self.__grid.shape[1] > y >= 0 and self.__grid[x][y] != '@':
            return True
        return False

    def execute_action(self, action):
        """
        :param action: 'north','south', 'east' or 'west' (chosen by user)
        :return: new current_state, reward
        """
        reward = self.get_reward(self.current_state, action)
        x, y = self.current_state
        actions = self.__get_actions(x, y)[action]
        probabilities = list(actions.values())
        possible_states = list(actions.keys())
        self.current_state = random.choices(possible_states, weights=probabilities)[0]
        return self.current_state, reward

    def is_episode_over(self):
        """
        :return whether the episode is over (reached end state)
        """
        if self.current_state == self.end_state:
            return True
        return False

    def new_episode(self):
        """
        :return whether the episode is over (reached end state)
        """
        self.current_state = self.start_state

    def create_q_values_table(self):
        """
        :return: dictionary in the form (x,y): q value of each direction
        """
        valid_locations = {}
        for i in range(self.__grid.shape[0]):
            for j in range(self.__grid.shape[1]):
                if self.__grid[i][j] != '@' and (i, j) != self.end_state:
                    valid_locations[(i, j)] = {'south': 0, 'north': 0, 'west': 0, 'east': 0}
        return valid_locations

    def create_policy_array(self):
        """
        :return: dictionary in the form (x,y): 'direction'
        """
        valid_locations = {}
        for i in range(self.__grid.shape[0]):
            for j in range(self.__grid.shape[1]):
                if self.__grid[i][j] != '@' and (i, j) != self.end_state:
                    valid_locations[(i, j)] = 'south'
        return valid_locations

    def pretty_print_policy(self, policy):
        pretty_policy = np.zeros(self.__grid.shape, 'U4')
        for state, direction in policy.items():
            if direction == 'east':
                pretty_policy[state[0]][state[1]] = '→'
            if direction == 'west':
                pretty_policy[state[0]][state[1]] = '←'
            if direction == 'north':
                pretty_policy[state[0]][state[1]] = '↑'
            if direction == 'south':
                pretty_policy[state[0]][state[1]] = '↓'
        for i, row in enumerate(pretty_policy):
            for j, cell in enumerate(row):
                if cell == '':
                    if self.end_state == (i, j):
                        pretty_policy[i][j] = 'G'
                    else:
                        pretty_policy[i][j] = '@'
        print(pd.DataFrame(pretty_policy, columns=[i for i in range(len(pretty_policy))]).to_markdown())
