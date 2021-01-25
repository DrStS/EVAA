import csv                  # csv file read
from math import isclose    # comparison

class CompareCsv:
    '''
        Compares the last row of two csv files
    '''

    def __init__(self):
        pass

    def _get_end_state_from_csv(self, path_to_csv_file):
        '''
            read in a .csv file and return the final row as dict
        '''
        with open(path_to_csv_file, 'r') as read_obj:
            dict_reader = csv.DictReader(read_obj, skipinitialspace=True)
            
            # all numeric values from the csv file will be of type string now
            list_of_dict = list(dict_reader)

            # cast the values from type string to type float
            list_of_dict = [dict([a, float(x)] for a, x in b.items()) for b in list_of_dict]

            # cast from OrderedDict to standard dict; ordering is detrimental for comparison
            return dict(list_of_dict[-1])

    def compare_end_state(self, path_to_ref_file, path_to_output_file, abstol, reltol):
        '''
            read in two csv files (really .csv files) and compare 
            the final row
        '''
        ref_end_state   = self._get_end_state_from_csv(path_to_ref_file)

        # end_state       = self._get_end_state_from_csv(path_to_output_file)
        # MOCK BEGIN
        end_state       = {'t': 1.000000, 'z_1': -0.006313108839, 'r_x': 0.003470642259, 'r_y': 0.001363759628, 'z_6': -0.002241651062, 'z_7': 0.000000000000, 'z_8': -0.000729903637, 'z_9': 0.000000000000, 'z_10': -0.000642214203, 'z_11': 0.000000000000, 'z_12': -0.002915203338, 'z_13': 0.000000000000}
        # MOCK END

        for key in ref_end_state.keys():
            if not(isclose(ref_end_state[key], end_state[key], abs_tol=float(abstol), rel_tol=float(reltol))):
                raise Exception('End states of ref file and output file are not equal for at least this key: ' + key + '(ref = ' + str(ref_end_state[key]) + ', output = ' + str(end_state[key]) + ')')