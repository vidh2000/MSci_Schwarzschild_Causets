import numpy as np
import cProfile
import pstats

def profiler(function,Nshow=None,args=None):
    """
    Profiles the function "function" with arguments
    "args" and returns the output of function(args)
    if args =! None.
    Processed are sorted by the TIME spent inside each of them
    and printed out. Shows "Nshow" operations
    """
    if args == None:
        with cProfile.Profile() as pr:    
            output = function()
            stats = pstats.Stats(pr)
            stats.sort_stats(pstats.SortKey.TIME)
            if Nshow == None:
                stats.print_stats()
            else:
                stats.print_stats(Nshow)    
    else:
        with cProfile.Profile() as pr:    
            output = function(args)
            stats = pstats.Stats(pr)
            stats.sort_stats(pstats.SortKey.TIME)
            if Nshow == None:
                stats.print_stats()
            else:
                stats.print_stats(Nshow)
    return output