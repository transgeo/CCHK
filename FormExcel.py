import pandas as pd
import numpy as np
import math
import copy
# ######################################################################################################################
# ######### Input data #########
station = pd.read_excel('HKMap.xlsx', sheet_name='Station_CC')
link = pd.read_excel('HKMap.xlsx', sheet_name='Link_CC')
trans = pd.read_excel('HKMap.xlsx', sheet_name='Transfer_station')
line_im = pd.read_excel('HKMap.xlsx', sheet_name='LineImportance')
LineColor = pd.read_excel('HKMap.xlsx', sheet_name='LineColor')
# ######### Setting parameters #########
line_num = 11
# Central point
central = np.array([114.180726, 22.334370])  # Kowloon Tong station
# Ring lines
ring_lines = np.array([1, 2, 3, 4, 5, 8, 10])
# Radial_lines
radial_lines = np.array([-1.0, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0.0,
                          0.125, 0.25, 0.375, 0.5, 0.675, 0.75, 0.875]) * math.pi
# ######################################################################################################################


# ######### Mapping to ring or radial lines ############################################################################
# ### MappingProcedure
# Output: MappingResult (array: station * (ID(started from 0) + ring line size + radial line size))
#         Indicating four possible mapping possible for each station
def MappingProcedure(r, theta, ring_lines, radial_lines):
    dr = None
    dtheta = None
    for i in range(ring_lines.size):
        if dr is None:
            dr = r - ring_lines[i]
        else:
            dr = np.vstack((dr, r - ring_lines[i]))
    for i in range(radial_lines.size):
        if dtheta is None:
            dtheta = theta - radial_lines[i]
        else:
            dtheta = np.vstack((dtheta, theta - radial_lines[i]))
    dr_t = copy.deepcopy(dr)
    dr_t[dr_t > 0] = 0
    dr_t[dr_t < 0] = 1
    dr_one = np.sum(dr_t, axis=0)
    dr_one[dr_one > 0] = 1  # dr_one: 0, located on the outside ring line
    dr_outer = copy.deepcopy(dr)
    dr_inter = copy.deepcopy(dr)
    dr_outer[dr_outer > 0] = -500000
    dr_inter[dr_inter < 0] = 500000
    dr_outer_loc = dr_outer.argmax(axis=0)  # 0, located on the outside ring line
    dr_inter_loc = dr_inter.argmin(axis=0)
    dtheta_clockwise = copy.deepcopy(dtheta)
    dtheta_counterclockwise = copy.deepcopy(dtheta)
    dtheta_clockwise[dtheta_clockwise < 0] = 4
    dtheta_counterclockwise[dtheta_counterclockwise > 0] = -4
    dtheta_clockwise_loc = dtheta_clockwise.argmin(axis=0)
    dtheta_counterclockwise_loc = dtheta_counterclockwise.argmax(axis=0)
    # ############# Adjusting
    station_ring = np.zeros((r.size, ring_lines.size), int)
    station_radial = np.zeros((r.size, radial_lines.size), int)
    for i in range(r.size):
        if dr_inter_loc[i] == (ring_lines.size - 1):
            station_ring[i][dr_inter_loc[i]] = 1
        elif r[i] == 0:  # origin node
            station_radial[i][:] = 1
        else:
            station_ring[i][dr_outer_loc[i]] = 1
            station_ring[i][dr_inter_loc[i]] = 1
            station_radial[i][dtheta_clockwise_loc[i]] = 1
            station_radial[i][dtheta_counterclockwise_loc[i]] = 1
    ID = np.linspace(0, r.size - 1, r.size, dtype=int).reshape([r.size, 1])
    MappingResult = np.hstack((ID, station_ring, station_radial))
    return MappingResult
# ######################################################################################################################

# ### TransferStationInformation  ######################################################################################
# Output: TM, list (bends information for each line for one transfer station)
#         TMM, (array, (# of transfer stations) * (transfer ID + line ID + station ID + ring / radial lines))
#         MP_update
def TransferStationInformation(MP, Line, transfer, line_im, ring_lines, radial_lines):
    TM = []
    for i in range(transfer.shape[0]):  # i for each transfer station
        ts = MP[transfer[i][2]-1, :]
        ts_in = np.where(ts[1:(ring_lines.size+radial_lines.size + 1)] == 1)
        TM_l = np.zeros((transfer[i][0] + 1, ts_in[0].size + 1), int)
        TM_l[0][1:ts_in[0].size + 1] = ts_in[0] + 1
        for j in range(transfer[i][0]):  # j for each line for transfer station i
            TM_l[j+1][0] = transfer[i][1+2*j]
            for t in range(ts_in[0].size):  # t for each possible mapping point
                MP_t = copy.deepcopy(MP)
                MP_t[transfer[i][2+2*j] - 1, 1:(ring_lines.size+radial_lines.size + 1)] = 0
                MP_t[transfer[i][2+2*j] - 1, ts_in[0][t]+1] = 1
                line_in = np.where(Line == transfer[i][1+2*j])
                line_index = line_in[0]
                LI = MP_t[line_index, :]
                LI_MP = AdjustingOneLine(LI, ring_lines, radial_lines)
                SI = SectionInformation(LI_MP)
                TM_l[j + 1][t+1] = SI.shape[0]
        TM.append(TM_l)
    TMM = np.zeros((np.sum(transfer[:, 0]), MP.shape[1]+3), int)
    m = 0
    for i in range(transfer.shape[0]):
        TT = TM[i]
        for j in range(transfer[i][0]):  # j for each line for transfer station i
            for k in range(TT.shape[1]-1):
                TMM[m][0] = i  # transfer ID
                TMM[m][1] = transfer[i][1 + 2 * j]  # line ID
                TMM[m][2] = transfer[i][2 + 2 * j]  # station ID
                TMM[m][TT[0][k+1]+2] = TT[j+1][k+1]  # station ID
            m = m + 1
    for i in range(line_num):
        LTL_in = np.where(TMM[:, 1] == line_im['Line ID'][i])
        LTL_index = LTL_in[0]
        LTL = TMM[LTL_index, :]
        LTL_t = copy.deepcopy(LTL)
        LTL_t[LTL_t == 0] = 10
        LTL_t[:, 3:(ring_lines.size + radial_lines.size + 4)] = LTL_t[:, 3:(ring_lines.size+radial_lines.size + 4)] + 1
        bb = np.argmin(LTL_t[:, 3:(ring_lines.size + radial_lines.size + 4)], axis=1)
        for j in range(bb.shape[0]):
            LTL_t[j, 3+bb[j]] = LTL_t[j, 3+bb[j]] - 1
        LTL_min = LTL_t[:, 3:(ring_lines.size+radial_lines.size + 4)].min(axis=1)
        LTL_t[:, 3:(ring_lines.size+radial_lines.size + 4)][LTL_t[:, 3:(ring_lines.size+radial_lines.size + 4)] > np.max(LTL_min)] = 0
        for j in range(LTL.shape[0]):
            TMM[TMM[:, 0] == LTL[j, 0], 3:(ring_lines.size+radial_lines.size + 4)] = LTL_t[j, 3:(ring_lines.size+radial_lines.size + 4)]
    MP_update = copy.deepcopy(MP)
    for i in range(TMM.shape[0]):
        a = TMM[i, 3:(ring_lines.size+radial_lines.size + 3)]
        a[a > 0] = 1
        MP_update[TMM[i][2]-1, 1:(ring_lines.size+radial_lines.size + 1)] = a
    return TM, TMM, MP_update
# ######################################################################################################################

# ######### Minimizing the bends for each line #########################################################################
# AdjustingEachLine  ---> MP
#   AdjustingOneLine  ---> LI (only one line information from MP)
#       SectionInformation  ---> SI (link section based on LI)
#       SelectingSection  ---> LI, LI_index, flag. LI:the update LI by selecting the link section with the most stations

# ### AdjustingEachLine
# Output: MappingResult (array: station * (ID(started from 0)+ ring line size + radial line size))
# Used function: SectionInformation, AdjustingOneLine
#       Finding the optimal location of each station
def AdjustingEachLine(Line, line_num, MP, ring_lines, radial_lines):
    MP_ae = copy.deepcopy(MP)
    for i in range(line_num):
        line_in = np.where(Line == i + 1)
        line_index = line_in[0]
        LI = MP[line_index, :]
        LI_MP = AdjustingOneLine(LI, ring_lines, radial_lines)
        MP_ae[line_index, :] = LI_MP
    return MP_ae
# ### AdjustingOneLine
# Output: LI MappingResult for one line
def AdjustingOneLine(LI, ring_lines, radial_lines):     # SI: section information for one line
    flag = 0
    LI_index = copy.deepcopy(LI)
    LI_ao = copy.deepcopy(LI)
    while flag == 0:
        SI = SectionInformation(LI_index)
        if SI.shape[0] == 0:
            break
        else:
            LI_ao, LI_index, flag = SelectingSection(LI_ao, LI_index, SI, ring_lines, radial_lines)
    return LI_ao
# ### SectionInformation
# Output: np.array(section_information)
#         (array: # of line sections (ring or radial line) * (ring/radial ID (column of MappingResult)
#                                                             + # of stations + start station ID + end station ID
#                                                       + start station ID (in one line) + end station ID (in one line))
# Input: MappingResult for one line
def SectionInformation(LI):     # LI: Line information for one line
    section_information = []
    for j in range(1, LI.shape[1]):
        index = 0
        for k in range(LI.shape[0]):
            if LI[k][j] == 1:
                if index == 0:
                    index = 1
                    start = LI[k][0]
                    start_l = k
                    end = LI[k][0]
                    end_l = k
                    n = 1
                else:
                    end = LI[k][0]
                    end_l = k
                    n = n + 1
            if LI[k][j] == 0 and index == 1:
                index = 0
                section_information.append([j, n, start, end, start_l, end_l])
            if index == 1 and k == LI.shape[0] - 1:
                index = 0
                section_information.append([j, n, start, end, start_l, end_l])
    return np.array(section_information)
# ### SelectingSection
def SelectingSection(LI, LI_index, SI, ring_lines, radial_lines):     # SI: section information for one line
    index = SI.argmax(axis=0)
    flag = 0
    LI_ss = copy.deepcopy(LI)
    LI_index_ss = copy.deepcopy(LI_index)
    if SI[index[1]][1] >= 1:
        LI_ss[SI[index[1]][4]:SI[index[1]][5]+1, 1:(ring_lines.size+radial_lines.size + 1)] = 0
        LI_ss[SI[index[1]][4]:SI[index[1]][5]+1, SI[index[1]][0]] = 1
    else:
        flag = 1
    LI_index_ss[SI[index[1]][4]:SI[index[1]][5] + 1, 1:(ring_lines.size+radial_lines.size + 1)] = 0
    return LI_ss, LI_index_ss, flag
# ######################################################################################################################

# ### LineInformation ##################################################################################################
# LineSection:  line section information for each line, cut by the transfer station
#               # of stations + ori_ID + des_ID + ori_trans + des_trans
# LinkSection:  link section information for each line section, link types between the transfer stations
#               type + radius / angle + theta 1 / r1 + theta 2 / r2 + # of stations
# LinkSectionA:  Adjusted one, the link section are continuous
#               type + radius / angle + theta 1 / r1 + theta 2 / r2 + # of stations
def LineInformation(Line, line_num, MP, r, theta, ring_lines, radial_lines):
    MP_up = np.where(MP[:, 1:(ring_lines.size + radial_lines.size + 1)] == 1)
    MP_sim = MP_up[1]
    StaLoc = np.zeros((MP.shape[0], 4), int)
    for i in range(MP.shape[0]):
        StaLoc[i][0] = i
        StaLoc[i][3] = MP_sim[i]
        if MP_sim[i] < ring_lines.size:
            StaLoc[i][1] = 1
            StaLoc[i][2] = MP_sim[i]
        else:
            StaLoc[i][1] = 2
            StaLoc[i][2] = MP_sim[i] - ring_lines.size

    LineSection = []
    for i in range(line_num):
        LS = []
        line_in = np.where(Line == i + 1)
        LI = StaLoc[line_in[0], :]
        for j in range(LI.shape[0]):
            if j == 0:
                link_ori = LI[j][0]
                link_des = LI[j][0]
                link_type = LI[j][1]
                ringorradial_ID = LI[j][2]
                ringandradial_ID = LI[j][3]
            elif j != LI.shape[0]-1:
                if LI[j][3] == ringandradial_ID:
                    link_des = LI[j][0]
                else:
                    if link_type == 1:
                        LS.append([link_type, link_ori, link_des,
                                   ring_lines[ringorradial_ID], theta[link_ori], theta[link_des]])
                    else:
                        LS.append([link_type, link_ori, link_des,
                                   radial_lines[ringorradial_ID], r[link_ori], r[link_des]])
                    link_ori = LI[j][0]
                    link_des = LI[j][0]
                    link_type = LI[j][1]
                    ringorradial_ID = LI[j][2]
                    ringandradial_ID = LI[j][3]
            elif j == LI.shape[0]-1:
                if LI[j][3] == ringandradial_ID:
                    link_des = LI[j][0]
                    if link_type == 1:
                        LS.append([link_type, link_ori, link_des,
                                   ring_lines[ringorradial_ID], theta[link_ori], theta[link_des]])
                    else:
                        LS.append([link_type, link_ori, link_des,
                                   radial_lines[ringorradial_ID], r[link_ori], r[link_des]])
                else:
                    if link_type == 1:
                        LS.append([link_type, link_ori, link_des,
                                   ring_lines[ringorradial_ID], theta[link_ori], theta[link_des]])
                    else:
                        LS.append([link_type, link_ori, link_des,
                                   radial_lines[ringorradial_ID], r[link_ori], r[link_des]])
                    if LI[j][1] == 1:
                        LS.append([LI[j][1], LI[j][0], LI[j][0],
                                   ring_lines[LI[j][2]], theta[LI[j][0]], theta[LI[j][0]]])
                    else:
                        LS.append([LI[j][1], LI[j][0], LI[j][0],
                                   radial_lines[LI[j][2]], r[LI[j][0]], r[LI[j][0]]])
        LineSection.append(LS)

    LineSectionA = []
    for i in range(line_num):
        LSA = []
        LS = LineSection[i]
        if LS.__len__() == 1:
            LineSectionA.append(LineSection[i])
        else:
            for j in range(LS.__len__()-1):
                if LS[j][0] == 1 and LS[j+1][0] == 1:
                    LSA.append([LS[j][0], LS[j][1], LS[j][2], LS[j][3],
                                LS[j][4], LS[j+1][4]])
                    LSA.append([2, LS[j][2], LS[j+1][1], LS[j+1][4],
                                LS[j][3], LS[j + 1][3]])
                if LS[j][0] == 1 and LS[j+1][0] == 2:
                    LSA.append([LS[j][0], LS[j][1], LS[j][2], LS[j][3],
                                LS[j][4], LS[j+1][3]])
                    LS[j + 1][4] = LS[j][3]
                if LS[j][0] == 2 and LS[j+1][0] == 1:
                    LSA.append([LS[j][0], LS[j][1], LS[j][2], LS[j][3],
                                LS[j][4], LS[j+1][3]])
                    LS[j+1][4] = LS[j][3]
                if LS[j][0] == 2 and LS[j + 1][0] == 2:
                    LSA.append([LS[j][0], LS[j][1], LS[j][2], LS[j][3],
                                LS[j][4], LS[j+1][4]])
                    LSA.append([1, LS[j][2], LS[j+1][1], LS[j+1][4],
                                LS[j][3], LS[j + 1][3]])
            LSA.append(LS[LS.__len__()-1])
            LineSectionA.append(LSA)
    return LineSectionA, StaLoc
# ######################################################################################################################

# ######################################################################################################################
# ######### Initial data processing #########
# Parameters for unit conversion
Ex = 100  # of Parameters for unit conversion of longitude and Cartesian coordinate
Ey = 100  # of Parameters for unit conversion of latitude and Cartesian coordinate
Ed = 1
# Initial data processing
Longitude = np.array(station['x'])
Latitude = np.array(station['y'])
Line = np.array(station['L'])
StationTrans = np.array(station['Trans'])
Stationcolor = station['color']
link_ori = np.array(link['Ori'])
link_des = np.array(link['Des'])
link_color = link['Color']
trans1 = np.array(trans)
transfer = trans1.astype(int)

# ######### Step 1 Determining the Cartesian coordinate #########
x = (Longitude - central[0]) * Ex
y = (Latitude - central[1]) * Ey

# ######### Step 2 Converting the Cartesian coordinate to Polar coordinate #########
r = np.sqrt(x ** 2 + y ** 2)
theta = np.arctan2(y, x)

# ######### Step 3 Mapping to ring line or radial line #########
MP = MappingProcedure(r, theta, ring_lines, radial_lines)
# ######### Step 4 determining the mapping point #########
# fixed transfer station
# minimize the number of bends
TI, TL, MP2 = TransferStationInformation(MP, Line, transfer, line_im, ring_lines, radial_lines)
MP3 = AdjustingEachLine(Line, line_num, MP2, ring_lines, radial_lines)
# output line information
LineS, StaLoc = LineInformation(Line, line_num, MP3, r, theta, ring_lines, radial_lines)

# write the station and link section information to excel file
writer = pd.ExcelWriter('Interactive_output_HK.xlsx')

station_output = pd.DataFrame(columns=['ID', 'lineID', 'Name', 'Trans', 'color',
                                       'r', 't', 'x', 'y', 'dx', 'dy'])

linesection_output = pd.DataFrame(columns=['ID', 'lineID',
                                           'Ori_name', 'Ori_ID', 'Des_name', 'Des_ID', 'type',
                                           'Ori_r', 'Ori_t', 'Ori_dr', 'Ori_dx ',
                                           'Des_r', 'Des_t', 'Des_dr', 'Des_dx',
                                           'Ori_x', 'Ori_y', 'Des_x', 'Des_y',
                                           'dt', 'LF', 'SF', 'color'])
# station info
r_out = []
t_out = []
for i in range(StaLoc.shape[0]):
    if StaLoc[i][1] == 1:
        r_out.append(ring_lines[StaLoc[i][2]])
        t_out.append(theta[i])
    else:
        t_out.append(radial_lines[StaLoc[i][2]])
        r_out.append(r[i])

for i in range(r_out.__len__()):
    b = pd.DataFrame([[i, Line[i], station['Station'][i],
                      station['Trans'][i], station['color'][i],
                      r_out[i], t_out[i],
                      r_out[i] * np.cos(t_out[i]) + 500,
                      500 - r_out[i] * np.sin(t_out[i]), 0, 0]],
                      columns=['ID', 'lineID', 'Name', 'Trans', 'color',
                               'r', 't', 'x', 'y', 'dx', 'dy'])
    station_output = station_output.append(b)

# line section info
ID = 1
for i in range(line_num):     # for each line
    Line_Section = LineS[i]
    Line_ID = i+1
    color = LineColor['Color'][i]
    for j in range(Line_Section.__len__()):   # for each section
        Ori_name = station['Station'][Line_Section[j][1]]
        Ori_ID = station['ID'][Line_Section[j][1]]
        Des_name = station['Station'][Line_Section[j][2]]
        Des_ID = station['ID'][Line_Section[j][2]]
        line_type = Line_Section[j][0]
        if line_type == 1:
            Ori_r = Line_Section[j][3]
            Ori_t = Line_Section[j][4]
            Des_r = Line_Section[j][3]
            Des_t = Line_Section[j][5]
        else:
            Ori_r = Line_Section[j][4]
            Ori_t = Line_Section[j][3]
            Des_r = Line_Section[j][5]
            Des_t = Line_Section[j][3]

        Ori_dr = 0
        Ori_dx = 0
        Des_dr = 0
        Des_dx = 0
        Ori_x = (Ori_r + Ori_dr) * np.cos(Ori_t) + Ori_dx + 500
        Ori_y = 500 - (Ori_r + Ori_dr) * np.sin(Ori_t)
        Des_x = (Des_r + Des_dr) * np.cos(Des_t) + Des_dx + 500
        Des_y = 500 - (Des_r + Des_dr) * np.sin(Des_t)
        dt = Des_t - Ori_t
        LF = 0
        if dt < 0:
            SF = 1
        else:
            SF = 0
        bb = pd.DataFrame([[ID, Line_ID, Ori_name, Ori_ID, Des_name, Des_ID, line_type,
                                   Ori_r, Ori_t, Ori_dr, Ori_dx,
                                   Des_r, Des_t, Des_dr ,Des_dx,
                                   Ori_x, Ori_y, Des_x, Des_y,
                                   dt, LF, SF,color]],
                                   columns=['ID', 'lineID',
                                            'Ori_name', 'Ori_ID', 'Des_name', 'Des_ID', 'type',
                                            'Ori_r', 'Ori_t', 'Ori_dr', 'Ori_dx ',
                                            'Des_r', 'Des_t', 'Des_dr', 'Des_dx',
                                            'Ori_x', 'Ori_y', 'Des_x', 'Des_y',
                                            'dt', 'LF', 'SF', 'color'])
        linesection_output = linesection_output.append(bb)
        ID = ID + 1


station_output.to_excel(writer, 'station')
linesection_output.to_excel(writer, 'linesection')
writer.save()
writer.close()
