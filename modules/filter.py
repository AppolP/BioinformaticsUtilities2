from typing import Union

def is_in_bounds(sequence: str, range: Union[tuple[int, int], int, float], *args: str) -> bool:
    '''
    Check if sequence length is in given bounds
    
    Arguments:
    sequence: str, sequence to process
    range: tuple, int or float, a range the sequence belongs to or a value the sequence is below 
    args: str, substrings to be counted and matched to the specified range
    '''
    if isinstance(range, (int, float)):
        low_bound, upper_bound = 0, range
    else:
        low_bound, upper_bound = range
    args_count = 0
    for arg in args:
        args_count += sequence.count(arg)  
    args_percentage = args_count/len(sequence)*100
    return low_bound <= args_percentage <= upper_bound

def is_qualified(seq_quality: str, quality_threshold: int) -> bool:
    '''
    Check if sequence quality under the given threshold
    
    Arguments: 
    seq_quality: str, quality sequence of the read sequence 
    quality_threshold: int, the value bordering acceptable quality value
    '''
    q_score = sum([ord(el)-33 for el in seq_quality])  
    return 10**(-q_score/10) >= quality_threshold