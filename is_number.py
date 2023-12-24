# This function check that our inputs for the filling percent of the vessel and the index for the hoppers are numbers, and are processable
def is_number(k,fill_percent):
    S = 0                             # If S remains 0, it menas that the information that user entered are processable. If it becomes 1, the inputs are either not ineger or not within the range.
    try:
        k =int(k)                   # The inputs are always considered as strings by Python. That is why we need to change k to a integer
    except ValueError:
        S = 1
    if (S==0 and (k > 6 or k < 1)):     # checking if the number is between 1 and 6
        S = 1
        #sys.exit(1)

    # checking if the input is a number
    try:
        fill_percent = float(fill_percent)  # The inputs are always considered as strings by Python. That is why we need to change fill_percent to a float
    except ValueError:
        S = 1
        #sys.exit(1)
    # checking if the number is between 0 and 100
    if (S==0 and (fill_percent > 100 or fill_percent < 0)):
        S = 1
        #sys.exit(1)

    return S, k, fill_percent
