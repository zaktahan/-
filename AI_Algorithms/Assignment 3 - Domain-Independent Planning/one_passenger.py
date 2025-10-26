from __future__ import print_function
from pyddl import Domain, Action, neg
from planner import planner

def create_domain_one_passenger():
    domain = Domain((
        Action(
            'move-up',  
            parameters=(
                ('taxi', 't'),
                ('position', 'px'),  # Current location on the x-axis
                ('position', 'py'),  # Current location on the y-axis
                ('position', 'by'),  # New location on the y-axis
            ),
            preconditions=(
                ('dec', 'py', 'by'),
                ('at', 't', 'px', 'py'),
            ),
            effects=(
                neg(('at', 't', 'px', 'py')),
                ('at', 't', 'px', 'by'),
            ),
        ),
        Action(
            'move-down',
            parameters=(
                ('taxi', 't'),
                ('position', 'px'),
                ('position', 'py'),
                ('position', 'by'),
            ),
            preconditions=(


            ),
            effects=(


            ),
        ),
        Action(
            'move-left',
            parameters=(
                ('taxi', 't'),
                ('position', 'px'),
                ('position', 'py'),
                ('position', 'bx'),
            ),
            preconditions=(


            ),
            effects=(


            ),
        ),
        Action(
            'move-right',
            parameters=(
                ('taxi', 't'),
                ('position', 'px'),
                ('position', 'py'),
                ('position', 'bx'),
            ),
            preconditions=(


            ),
            effects=(


            ),
        ),
        Action(
            'pick-up',
            parameters=(
                ('taxi', 't'),
                ('position', 'px'),
                ('position', 'py'),
                ('passenger', 'p'),
            ),
            preconditions=(


            ),
            effects=(


            ),
        ),
        Action(
            'put-down',
            parameters=(
                ('taxi', 't'),
                ('position', 'px'),
                ('position', 'py'),
                ('passenger', 'p'),
            ),
            preconditions=(


            ),
            effects=(


            ),
        ),
    ))
    return domain

