import numpy as np
import cProfile
import pstats

def profiler(function, Nshow=None, *args, **kwargs):
    """
    Profiles the function "function" with arguments
    "*args" and key arguments "**kwargs". Prints stats, sorted
    by the TIME spent inside each.\n
    
    Parameters
    ----------
    - function : Callable.
    
    - Nshow : int\n
        Print"Nshow" operations if specified. Otherwise, if Nshow=None
        (Default) print all. 
    
    - *args, **kwargs of function

    Returns
    --------
    output: whatever
        The output of function(*args, **kwargs)
    """
    with cProfile.Profile() as pr:    
        output = function(*args, **kwargs)
        stats = pstats.Stats(pr)
        stats.sort_stats(pstats.SortKey.TIME)
        if Nshow == None:
            stats.print_stats()
        else:
            stats.print_stats(Nshow)  
    # if args == None:
    #     with cProfile.Profile() as pr:    
    #         output = function()
    #         stats = pstats.Stats(pr)
    #         stats.sort_stats(pstats.SortKey.TIME)
    #         if Nshow == None:
    #             stats.print_stats()
    #         else:
    #             stats.print_stats(Nshow)    
    # else:
    #     with cProfile.Profile() as pr:    
    #         output = function(args)
    #         stats = pstats.Stats(pr)
    #         stats.sort_stats(pstats.SortKey.TIME)
    #         if Nshow == None:
    #             stats.print_stats()
    #         else:
    #             stats.print_stats(Nshow)
    return output