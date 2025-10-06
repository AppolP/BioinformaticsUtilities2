def is_in_bounds(sequence: str, range: Union[tuple[int, int], int, float], *args: str) -> bool:
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
    q_score = sum([ord(el)-33 for el in seq_quality])  
    return 10**(-q_score/10) >= quality_threshold