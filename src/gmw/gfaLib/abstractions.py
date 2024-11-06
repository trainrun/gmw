from enum import Enum

class Orientation(Enum):

    FORWARD = '+'
    REVERSE = '-'
    ANY = '?'
    BOTH = '='


def reverse(orientation: Orientation) -> Orientation:

    return {
        Orientation.FORWARD: Orientation.REVERSE,
        Orientation.REVERSE: Orientation.FORWARD,
        Orientation.ANY: Orientation.ANY,
        Orientation.BOTH: Orientation.BOTH,
    }[orientation]


class GFAFormat(Enum):
    RGFA = 'rGFA'
    GFA1 = 'GFA1'
    GFA1_1 = 'GFA1.1'
    GFA1_2 = 'GFA1.2'
    GFA2 = 'GFA2'
    ANY = 'unknown'


class GFALine(Enum):
    SEGMENT = 'S'
    LINK = 'L'
    WALK = 'W'
    PATH = 'P'
    HEADER = 'H'
    COMMENT = '#'
    ANY = '?'
