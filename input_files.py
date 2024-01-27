
def input_files():

    import pandas as pd

    [df_FFC, df_wall] = input_files()
    # reading the sheet that has the stress walk and bulk density data vs major principal stress
    df_FFC = pd.read_excel('powder_rheometry_data2.xls', sheet_name='stress_walk')
    # reading the sheet that has the wall friction data (shear stress as a function of normal stress)
    df_wall = pd.read_excel('powder_rheometry_data2.xls', sheet_name='wall_friction_data')

    return df_FFC, df_wall