
import pandas as pd
import numpy as np

# read function
def read_temperature_data(file_path, column_name):
    df = pd.read_csv(file_path)
    df['time'] = pd.to_datetime(df['time'])
    df[column_name] = df[column_name].str.replace(" °C", "").astype(float)
    df["time_15min"] = df["time"].dt.floor("15min")
    return df.groupby("time_15min")[[column_name]].mean().reset_index()[column_name].values

# read temperature
Tem_ind = read_temperature_data("data/TemperaturesInd.csv", "Room 1 - Actual")
T_amb = read_temperature_data("data/TemperaturesOut.csv", "Outside")

# read network
network = pd.read_csv("data/network.csv")
connect1 = pd.read_excel("data/user_connect.xlsx")

# read price
price_data = pd.read_csv("data/price_data.csv")
filtered_data = price_data.query("year == 2024 & month == 2 & day == 1")

# read load
def read_power_data(file_path):
    df = pd.read_csv(file_path, index_col=0)
    df['time'] = pd.to_datetime(df['time'], format='%m/%d/%Y %I:%M %p')
    return df[df['time'].dt.date == pd.to_datetime("2024-02-01").date()].reset_index(drop=True)

LoadPower = read_power_data("data/UserPower.csv")
LoadReact = read_power_data("data/UserReactivePower.csv")

# pv factor
weather = pd.read_csv('data/weather.csv')
weather['timestamp'] = pd.to_datetime(weather['timestamp'])
solar = weather.loc[weather['timestamp'].dt.date == pd.Timestamp('2024-02-01').date(), ['timestamp', 'P_solar']]

# 15 min time index
time_index_15min = pd.date_range(start=solar['timestamp'].min(), 
                                 end=solar['timestamp'].max(), 
                                 freq='15min')

# reindex and interpolate
solar_15min = solar.set_index('timestamp').reindex(time_index_15min).interpolate()
solar_15min.reset_index(inplace=True)
solar_15min.rename(columns={'index': 'timestamp'}, inplace=True)
solar_output = solar_15min['P_solar'].values /1000


pv_factor = np.repeat([
    0, 0, 0, 0, 0, 0, 0, 0, 0.006, 0.053, 0.129, 0.179,
    0.166, 0.14, 0.094, 0.046, 0.007, 0, 0, 0, 0, 0, 0, 0
], 4)
pv_factor = pd.DataFrame(pv_factor, columns=["Quarter_Hourly_Data"])

# congestion set
congestion_limit = np.ones(96) * (-50 * 1e-3)
congestion_limit[40:48] = -100 * 1e-3  # 设定拥塞区间


class ModelInf:
    def __init__(self):
        # network parameter
        self.s_trafo = 0.23 * 0.16 * 3
        self.pf_pv_limit = 0.95
        self.tan_phi_pv = (np.sqrt(1 - self.pf_pv_limit**2)) / self.pf_pv_limit
        self.pf = 0.9
        self.tan_phi_load = (np.sqrt(1 - self.pf**2)) / self.pf
        self.v_ref = 0.23
        self.v_lb = 0.96
        self.v_ub = 1.04

        # HHP
        self.p_hp_max = 5e-3  # 5kw
        self.p_hp_min = 3e-3  # 1.8kw
        self.p_boil_max = 15e-3  # 15kw
        self.p_boil_min = 1e-3  # 1kw
        self.gas_LHV = 10.16  # kWh/m3
        self.COP = 4.5
        self.hp_own=connect1["HP"].tolist()

        # indoor thermal
        self.C_house = 10 + 0.2 + 30
        self.R_house = 1 / (1/7 + 1/6 + 1/2)
        self.x_vals = list(range(19, 24))
        self.y_vals = [12.3, 9.9, 8.2, 7.2, 6.6]
        self.Pn = len(self.x_vals)

        # price
        self.PPD_Price = 0.27
        self.ele_price = filtered_data['energy_price_full'].values
        self.gas_price = filtered_data['gas_price_full'].values

        # load data
        self.LoadPower = LoadPower
        self.LoadReact = LoadReact
        self.pvFactor = pv_factor
        self.solar_output = solar_output

        # congestiion
        self.congestion_limit = congestion_limit
        self.StartRe = 40
        self.EndRe = 48

        # network input
        self.network = network
        self.connect1 = connect1

        # temperature data
        self.T_amb = T_amb
        self.Tem_ind = Tem_ind




def get_model_inf():
    return ModelInf()




