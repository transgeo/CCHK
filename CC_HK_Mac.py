import pandas as pd
import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import seaborn as sns


station = pd.read_excel('HKMap.xlsx', sheet_name='Station_CC')
link = pd.read_excel('HKMap.xlsx', sheet_name='Link_CC')
trans = pd.read_excel('HKMap.xlsx', sheet_name='Transfer_station')
line_im = pd.read_excel('HKMap.xlsx', sheet_name='LineImportance')
line_color = pd.read_excel('HKMap.xlsx', sheet_name='LineColor')

# ######### Input data #########

# ######### Setting parameters #########
line_num = 11
# Central point
central = np.array([114.180726, 22.334370])  #
# Ring lines
ring_lines = np.array([1, 2, 3, 4, 5, 8, 10])
# 2, 5, 10, 15, 20, 25
# Radial_lines
radial_lines = np.array([-1.0, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0.0,
                          0.125, 0.25, 0.375, 0.5, 0.675, 0.75, 0.875]) * math.pi
# central = [114.186686, 22.299936]   # Hong Hom
# central = [114.194626, 22.302199]   # Whampoa
# central = [114.164061, 22.278969]   # Central
# central = [114.180726, 22.334370]   # Kowloon
# ring_lines = [1, 3, 6, 9]
# ring_lines = [1, 2, 3, 4, 5, 8, 10]
# radial_lines = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75]
# radial_lines = [-1, -0.833, -0.667, -0.5, -0.333, -0.167, 0, 0.167, 0.333, 0.5, 0.667, 0.833]
# radial_lines = [-1.0, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0.0,
#                          0.125, 0.25, 0.375, 0.5, 0.675, 0.75, 0.875]
# ######### Initial data processing #########


# Parameters for unit conversion
Ex = 100  # of Parameters for unit conversion of longitude and Cartesian coordinate
Ey = 100  # of Parameters for unit conversion of latitude and Cartesian coordinate
Ed = 1
# ######### Initial data processing #########
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
    dr_outer[dr_outer > 0] = -25
    dr_inter[dr_inter < 0] = 25
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

# ######### Step 3 Mapping to ring line or radial line #########
MP = MappingProcedure(r, theta, ring_lines, radial_lines)

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
# ### TransferStationInformation
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
                MP_t[transfer[i][2+2*j] - 1, 1:24] = 0
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
        LTL_min = LTL_t[:, 3:(ring_lines.size+radial_lines.size + 4)].min(axis=1)
        LTL[:, 3:(ring_lines.size+radial_lines.size + 4)][LTL[:, 3:(ring_lines.size+radial_lines.size + 4)] > np.max(LTL_min)] = 0
        for j in range(LTL.shape[0]):
            TMM[TMM[:, 0] == LTL[j, 0], 3:(ring_lines.size+radial_lines.size + 4)] = LTL[j, 3:(ring_lines.size+radial_lines.size + 4)]
    MP_update = copy.deepcopy(MP)
    for i in range(TMM.shape[0]):
        a = TMM[i, 3:(ring_lines.size+radial_lines.size + 3)]
        a[a > 0] = 1
        MP_update[TMM[i][2]-1, 1:(ring_lines.size+radial_lines.size + 1)] = a
    return TM, TMM, MP_update

# ######### Step 4 Adjusting #########
# ######### Step 4.1 Minimizing the bends #########
MP1 = AdjustingEachLine(Line, line_num, MP, ring_lines, radial_lines)
TI, TL, MP2 = TransferStationInformation(MP, Line, transfer, line_im, ring_lines, radial_lines)
MP3 = AdjustingEachLine(Line, line_num, MP2, ring_lines, radial_lines)
# ######### Step 4.2 Averaging the link length #########

# ### LineInformation
# LineSection:  line section information for each line, cut by the transfer station
#               # of stations + ori_ID + des_ID + ori_trans + des_trans
# LinkSection:  link section information for each line section, link types between the transfer stations
#               type + radius / angle + theta 1 / r1 + theta 2 / r2 + # of stations
# LinkSectionA:  Adjusted one, the link section are continuous
#               type + radius / angle + theta 1 / r1 + theta 2 / r2 + # of stations
def LineInformation(Line, line_num, MP, r, theta, ST, ring_lines, radial_lines):
    LineSection = []
    for i in range(line_num):
        LS = []
        line_in = np.where(Line == i + 1)
        line_index = line_in[0]
        TS_information = ST[line_index]
        NoS = 0
        for j in range(TS_information.size):
            if j == 0:
                ori_ID = line_index[j]
                des_ID = line_index[j]
                ori_trans = TS_information[j]
                des_trans = TS_information[j]
                NoS = 1
            elif j == TS_information.size - 1:
                des_ID = line_index[j]
                des_trans = TS_information[j]
                NoS = 1 + NoS
                LS.append([NoS, ori_ID, des_ID, ori_trans, des_trans])
            else:
                if TS_information[j] == 1:
                    des_ID = line_index[j]
                    des_trans = TS_information[j]
                    NoS = 1 + NoS
                    LS.append([NoS, ori_ID, des_ID, ori_trans, des_trans])
                    ori_ID = line_index[j]
                    ori_trans = TS_information[j]
                    NoS = 1
                else:
                    NoS = 1 + NoS
        LineSection.append(LS)
    LinkSectionIn = []
    for i in range(line_num):
        LS = LineSection[i]
        LinkS = []
        for j in range(LS.__len__()):
            LI = MP[LS[j][1]:LS[j][2]+1, :]
            LI_index = np.where(LI[:, 1:(ring_lines.size+radial_lines.size + 1)] == 1)[1]
            LinkIn = []
            for k in range(LS[j][0]):
                if k == 0:
                    ori = LS[j][1]
                    des = LS[j][1]
                    rl = LI_index[0]
                elif k == LS[j][0] - 1:
                    if LI_index[k] == rl:
                        des = LS[j][1] + k
                        LinkIn.append([rl, ori, des])
                    else:
                        LinkIn.append([rl, ori, des])
                        LinkIn.append([LI_index[k], LS[j][1] + k, LS[j][1] + k])
                else:
                    if LI_index[k] == rl:
                        des = LS[j][1] + k
                    else:
                        LinkIn.append([rl, ori, des])
                        ori = LS[j][1] + k
                        des = LS[j][1] + k
                        rl = LI_index[k]
            LinkS.append(LinkIn)
        LinkSectionIn.append(LinkS)
    LinkSection = []
    LinkSectionA = []
    for i in range(line_num):
        LS = LineSection[i]
        LSI = LinkSectionIn[i]
        link_sect = []
        link_sectA = []
        for j in range(LS.__len__()):
            LinkSec = []
            LSI_e = LSI[j]
            for k in range(LSI[j].__len__()):
                if LSI_e[k][0] < ring_lines.size:
                    LinkSec.append([1, ring_lines[LSI_e[k][0]], theta[LSI_e[k][1]], theta[LSI_e[k][2]], LSI_e[k][2] - LSI_e[k][1] + 1])
                    # type    r   theta 1    theta 2    number of stations
                else:
                    # type    theta    r1    r2    number of stations
                    LinkSec.append([2, radial_lines[LSI_e[k][0] - ring_lines.size], r[LSI_e[k][1]], r[LSI_e[k][2]], LSI_e[k][2] - LSI_e[k][1] + 1])
            LinkSecA = []
            if LSI[j].__len__() > 1:                    # adjusting
                for k in range(LSI[j].__len__()-1):
                    if LinkSec[k][0] == 1:
                        if LinkSec[k+1][0] == 1:
                            LinkSec[k][3] = LinkSec[k + 1][2]
                            LinkSecA.append(LinkSec[k])
                            LinkSecA.append([2, LinkSec[k][3], LinkSec[k][1], LinkSec[k+1][1], 0])
                        else:
                            if LinkSec[k+1][2] != 0:
                                LinkSec[k][3] = LinkSec[k + 1][1]
                                LinkSecA.append(LinkSec[k])
                                LinkSec[k+1][2] = LinkSec[k][1]
                            else:
                                LinkSec[k + 1][1] = LinkSec[k][3]
                                LinkSecA.append(LinkSec[k])
                                LinkSec[k + 1][2] = LinkSec[k][1]
                    else:
                        if LinkSec[k][2] != 0:
                            if LinkSec[k+1][0] == 1:
                                LinkSec[k][3] = LinkSec[k + 1][1]
                                LinkSecA.append(LinkSec[k])
                                LinkSec[k + 1][2] = LinkSec[k][1]
                            else:
                                if abs(LinkSec[k][1] - LinkSec[k+1][1]) - np.pi > 0.01:
                                    LinkSec[k][3] = LinkSec[k + 1][2]
                                    LinkSecA.append(LinkSec[k])
                                    LinkSecA.append([1, LinkSec[k][3], LinkSec[k][1], LinkSec[k + 1][1], 0])
                                else:
                                    LinkSec[k][3] = 0
                                    LinkSecA.append(LinkSec[k])
                                    LinkSec[k+1][2] = 0
                        else:
                            if k > 0:
                                if LinkSec[k - 1][0] == 1:
                                    LinkSec[k][1] = LinkSec[k - 1][3]
                                    LinkSec[k][2] = LinkSec[k - 1][1]
                                else:
                                    LinkSec[k][1] = LinkSec[k - 1][1]
                                    LinkSec[k][2] = LinkSec[k - 1][3]
                                if LinkSec[k + 1][0] == 1:
                                    LinkSec[k][3] = LinkSec[k + 1][1]
                                    LinkSecA.append(LinkSec[k])
                                    LinkSec[k + 1][2] = LinkSec[k][1]
                                else:
                                    if abs(LinkSec[k][1] - LinkSec[k + 1][1]) - np.pi < 0.01:
                                        LinkSec[k][3] = LinkSec[k + 1][2]
                                        LinkSecA.append(LinkSec[k])
                                        LinkSecA.append([1, LinkSec[k][3], LinkSec[k][1], LinkSec[k + 1][1], 0])
                                    else:
                                        LinkSec[k][3] = 0
                                        LinkSecA.append(LinkSec[k])
                                        LinkSec[k + 1][2] = 0
                            else:
                                if LinkSec[k + 1][0] == 1:
                                    LinkSec[k][3] = LinkSec[k + 1][1]
                                    LinkSec[k][1] = LinkSec[k + 1][2]
                                    LinkSecA.append(LinkSec[k])
                                    LinkSec[k + 1][2] = LinkSec[k][1]
                                else:
                                    LinkSec[k][3] = LinkSec[k + 1][2]
                                    LinkSec[k][1] = LinkSec[k + 1][1]
                                    LinkSecA.append(LinkSec[k])
                                    LinkSecA.append([1, LinkSec[k][3], LinkSec[k][1], LinkSec[k + 1][1], 0])
                LinkSecA.append(LinkSec[LSI[j].__len__()-1])
            else:
                LinkSecA = copy.deepcopy(LinkSec)
            link_sect.append(LinkSec)
            link_sectA.append(LinkSecA)
        LinkSection.append(link_sect)
        LinkSectionA.append(link_sectA)
    LinkSectionAL = []
    for i in range(line_num):
        LSA = copy.deepcopy(LinkSectionA[i])
        Lin = copy.deepcopy(LineSection[i])
        LSAal = []
        for j in range(LSA.__len__()):
            LSAL = LSA[j]
            LinL = Lin[j]
            sl = 0
            for k in range(LSAL.__len__()):
                if LSAL[k][0] == 1:
                    if LSAL[k][2] * LSAL[k][3]>0:
                        sl = sl + np.abs(LSAL[k][2] - LSAL[k][3]) * LSAL[k][1]
                    else:
                        if np.abs(LSAL[k][3]) > np.pi * 0.5:
                            dt = 2 * np.pi - np.abs(LSAL[k][3]) - np.abs(LSAL[k][2])
                        else:
                            dt = np.abs(LSAL[k][3]) + np.abs(LSAL[k][2])
                        sl = sl + dt * LSAL[k][1]
                else:
                    sl = sl + np.abs(LSAL[k][3] - LSAL[k][2])
            ds = sl / (LinL[0] - 1)
            if ds < 0.2:
                esl = 0.5 * (LinL[0] - 1) - sl
                if LinL[3] == 0:
                    if LSAL[0][0] == 1:
                        LSAL[0][2] = LSAL[0][2] - esl / LSAL[0][1] * (
                                    np.abs(LSAL[0][2] - LSAL[0][3]) / (LSAL[0][2] - LSAL[0][3]))
                    else:
                        LSAL[0][2] = LSAL[0][2] - esl * (
                                np.abs(LSAL[0][3] - LSAL[0][2]) / (LSAL[0][3] - LSAL[0][2]))
                if LinL[4] == 0:
                    if LSAL[LSAL.__len__()-1][0] == 1:
                        LSAL[LSAL.__len__() - 1][3] = LSAL[LSAL.__len__() - 1][3] + esl / LSAL[LSAL.__len__() - 1][
                            1] * (
                                                              np.abs(LSAL[LSAL.__len__() - 1][3] -
                                                                     LSAL[LSAL.__len__() - 1][2]) / (
                                                                      LSAL[LSAL.__len__() - 1][3] -
                                                                      LSAL[LSAL.__len__() - 1][2]))
                    else:
                        LSAL[LSAL.__len__() - 1][3] = LSAL[LSAL.__len__() - 1][3] + esl * (
                                np.abs(LSAL[LSAL.__len__() - 1][3] -
                                       LSAL[LSAL.__len__() - 1][2]) / (
                                        LSAL[LSAL.__len__() - 1][3] -
                                        LSAL[LSAL.__len__() - 1][2]))
            LSAal.append(LSAL)
        LinkSectionAL.append(LSAal)
    return LineSection, LinkSectionIn, LinkSection, LinkSectionA, LinkSectionAL

LineS, LinkSIn, LinkS, LinkA, LinkSA = LineInformation(Line, line_num, MP3, r, theta, StationTrans, ring_lines, radial_lines)

# ##### averagelength
# station: type, r, theta
def averagelength(LineSection, LinkSection):
    station = []
    LengthInfo = []
    for i in range(LineSection.__len__()):   # for each line
        LeS = LineSection[i]
        LkS = LinkSection[i]
        stationl = []
        LengthInfol = []
        for j in range(LeS.__len__()):      # for each section between the transfer stations
            tl = 0
            nl = LeS[j][0]
            LkSp = LkS[j]
            sl = np.zeros((1, LkSp.__len__()))
            stationll = []
            for k in range(LkSp.__len__()):   # tl, sl,   # line sections between transfers stations
                if LkSp[k][0] == 1:
                    if LkSp[k][3] * LkSp[k][2] > 0:
                        sl[0][k] = sl[0][k] + np.abs((LkSp[k][3] - LkSp[k][2]) * LkSp[k][1])
                    else:
                        if np.abs(LkSp[k][3]-LkSp[k][2]) > np.pi:
                            dt = 2 * np.pi - np.abs(LkSp[k][3]) - np.abs(LkSp[k][2])
                        else:
                            dt = np.abs(LkSp[k][3]-LkSp[k][2])
                        sl[0][k] = sl[0][k] + dt * LkSp[k][1]
                else:
                    sl[0][k] = sl[0][k] + np.abs(LkSp[k][3] - LkSp[k][2])
                tl = tl + sl[0][k]
            ds = tl / (nl - 1)
            bl = 0
            LengthInfol.append([tl, ds, nl])
            if LkSp[0][0] == 1:
                stationll.append([LkSp[0][0], LkSp[0][1], LkSp[0][2]])
            else:
                stationll.append([LkSp[0][0], LkSp[0][2], LkSp[0][1]])
            for k in range(LkSp.__len__()):
                rl = sl[0][k] + bl
                if LkSp[k][0] == 1:
                    ddt = ds / LkSp[k][1]
                    sst = bl / LkSp[k][1]
                    stat = LkSp[k][2]
                    endt = LkSp[k][3]
                    if (stat * endt < 0):
                        if stat < -0.5 * np.pi:
                            stat = stat + 2 * np.pi
                        if endt < -0.5 * np.pi:
                            endt = endt + 2 * np.pi
                    fi = np.abs(endt - stat)/(endt - stat)
                    while (rl - ds) > -0.001:
                        in_tt = stat+ddt*fi-sst*fi
                        if (in_tt < - np.pi):
                            in_tt = in_tt +2 * np.pi
                        if (in_tt > np.pi):
                            in_tt = in_tt - 2 * np.pi
                        stationll.append([LkSp[k][0], LkSp[k][1], in_tt])
                        rl = rl - ds
                        stat = stat+ddt*fi
                    bl = rl
                else:
                    star = LkSp[k][2]
                    fi = np.abs(LkSp[k][3] - LkSp[k][2]) / (LkSp[k][3] - LkSp[k][2])
                    while (rl - ds) > -0.001:
                        stationll.append([LkSp[k][0], star + ds * fi - bl * fi, LkSp[k][1]])
                        rl = rl - ds
                        star = star + ds * fi
                    bl = rl
            stationl.append(stationll)
        station.append(stationl)
        LengthInfo.append(LengthInfol)
    return station, LengthInfo

st, LenInfo = averagelength(LineS, LinkSA)
# ##### Overlap
# OV: type, # of lines used (current ring line), # of lines used (current radial line)
def Overlap(station, radial_lines, ring_lines):
    OV = []
    ring_used = np.zeros((ring_lines.size, 8))
    radial_used = np.zeros((1, radial_lines.size))
    for i in range(station.__len__()):
        stl = station[i]
        OVl = []
        ring_usedl = np.zeros((ring_lines.size, 8))
        radial_usedl = np.zeros((1, radial_lines.size))
        for j in range(stl.__len__()):
            stll = stl[j]
            OVll = []
            for k in range(stll.__len__()):
                if stll[k][0] == 1:
                    if sum(np.where(abs(ring_lines - stll[k][1]) < 0.001, 1, 0))==0:
                        OVll.append([stll[k][0], 0, 0])
                    else:
                        ring_index = np.where(abs(ring_lines - stll[k][1]) < 0.001)[0][0]
                        ring_section = int(stll[k][2] / np.pi * 4) + 4
                        ring_usedl[ring_index, ring_section] = 1
                        if ring_used[ring_index, ring_section] != 0:
                            OVll.append([stll[k][0], ring_used[ring_index, ring_section], 0])
                        else:
                            OVll.append([stll[k][0], 0, 0])
                else:
                    if sum(np.where(abs(radial_lines - stll[k][2])<0.001, 1, 0))==0:
                        OVll.append([stll[k][0], 0, 0])
                    else:
                        radial_index = np.where(abs(radial_lines - stll[k][2])<0.001)[0][0]
                        radial_usedl[0, radial_index] = 1
                        if radial_used[0, radial_index] != 0:
                            OVll.append([stll[k][0], 0, radial_used[0, radial_index]])
                        else:
                            OVll.append([stll[k][0], 0, 0])
            OVl.append(OVll)
        ring_used = ring_used + ring_usedl
        radial_used = radial_used + radial_usedl
        OV.append(OVl)
    return OV

OV = Overlap(st, radial_lines, ring_lines)

# ######################### Plot the map ###############
rr = []
tt = []
ty = []
rr_ov = []
tt_ov = []
for i in range(st.__len__()):
    for j in range(st[i].__len__()):
        if j == 0:
            for k in range(st[i][j].__len__()):
                rr.append(st[i][j][k][1])
                tt.append(st[i][j][k][2])
                rr_ov.append(OV[i][j][k][1])
                tt_ov.append(OV[i][j][k][2])
                ty.append(st[i][j][k][0])
        else:
            for k in range(1, st[i][j].__len__()):
                rr.append(st[i][j][k][1])
                tt.append(st[i][j][k][2])
                rr_ov.append(OV[i][j][k][1])
                tt_ov.append(OV[i][j][k][2])
                ty.append(st[i][j][k][0])

x = []
y = []
dr = 0.2
dt = -0.2
for i in range(rr.__len__()):
    if rr_ov[i] != 0:
        r_a = rr[i] + rr_ov[i] * dr
        x.append(r_a * np.cos(tt[i]))
        y.append(r_a * np.sin(tt[i]))
    elif tt_ov[i] != 0:
        x.append(rr[i] * np.cos(tt[i]))
        y.append(rr[i] * np.sin(tt[i]) + tt_ov[i]*dt)
    else:
        x.append(rr[i] * np.cos(tt[i]))
        y.append(rr[i] * np.sin(tt[i]))

plt.scatter(x, y, c='b', s=5)

for i in range(link_ori.size):
    ori = link_ori[i]-1
    des = link_des[i]-1
    co = link_color[i]
    if rr[ori] == rr[des]:
        if rr_ov[ori] != 0:
            r_a = rr[ori] + rr_ov[ori] * dr
        else:
            r_a = rr[ori]
        if tt[ori]*tt[des] < 0:
            if tt[ori] < -0.5*np.pi:
                tt[ori] = tt[ori] + 2*np.pi
            if tt[des] < -0.5 * np.pi:
                tt[des] = tt[des] + 2 * np.pi
        t = np.linspace(tt[ori], tt[des], 100)
        x1, y1 = np.cos(t) * r_a, np.sin(t) * r_a
        plt.plot(x1, y1, color=co)
    elif tt[ori] == tt[des]:
        if tt_ov[ori] != 0:
            x1 = rr[ori] * np.cos(tt[ori]), rr[des] * np.cos(tt[des])
            y1 = rr[ori] * np.sin(tt[ori]) + tt_ov[ori]*dt, rr[des] * np.sin(tt[des]) + tt_ov[des]*dt
            plt.plot(x1, y1, color=co)
        else:
            x1 = rr[ori] * np.cos(tt[ori]), rr[des] * np.cos(tt[des])
            y1 = rr[ori] * np.sin(tt[ori]), rr[des] * np.sin(tt[des])
            plt.plot(x1, y1, color=co)
    else:
        if rr_ov[ori] != 0:
            ori_r = rr[ori] + rr_ov[ori] * dr
        else:
            ori_r = rr[ori]
        if rr_ov[des] != 0:
            des_r = rr[des] + rr_ov[des] * dr
        else:
            des_r = rr[des]
        ori_t = tt[ori]
        des_t = tt[des]
        if ty[ori] == 1:
            r_in = ori_r
            t_in = des_t
            if ori_t*t_in<0:
               if ori_t < -0.5 * np.pi:
                   ori_t = ori_t + 2 * np.pi
               if t_in < -0.5 * np.pi:
                   t_in = t_in + 2 * np.pi
            t = np.linspace(ori_t, t_in, 100)
            x1, y1 = np.cos(t) * ori_r, np.sin(t) * r_in
            plt.plot(x1, y1, color=co)
            if tt_ov[des] != 0:
                x1 = r_in * np.cos(t_in), des_r * np.cos(des_t)
                y1 = r_in * np.sin(t_in) + tt_ov[des] * dt, des_r * np.sin(des_t)+ tt_ov[des] * dt
                y2 = r_in * np.sin(t_in), r_in * np.sin(t_in) + tt_ov[ori] * dt
                x2 = des_r * np.cos(des_t), des_r * np.cos(des_t)
                plt.plot(x2, y2, color=co)
            else:
                x1 = r_in * np.cos(t_in), des_r * np.cos(des_t)
                y1 = r_in * np.sin(t_in), des_r * np.sin(des_t)
            plt.plot(x1, y1, color=co)
        else:
            if ori_r !=0 and des_r != 0:
                r_in = des_r
                t_in = ori_t
                if tt_ov[ori] != 0:
                    x1 = ori_r * np.cos(ori_t), r_in * np.cos(t_in)
                    y1 = ori_r * np.sin(ori_t) + tt_ov[ori] * dt, r_in * np.sin(t_in) + tt_ov[ori] * dt
                    x2 = r_in * np.cos(t_in), r_in * np.cos(t_in)
                    y2 = r_in * np.sin(t_in) + tt_ov[ori] * dt, r_in * np.sin(t_in)
                    plt.plot(x2, y2, color=co)
                else:
                    x1 = ori_r * np.cos(ori_t), r_in * np.cos(t_in)
                    y1 = ori_r * np.sin(ori_t), r_in * np.sin(t_in)
                plt.plot(x1, y1, color=co)
                if des_t * t_in < 0:
                    if des_t < -0.5 * np.pi:
                        des_t = des_t + 2 * np.pi
                    if t_in < -0.5 * np.pi:
                        t_in = t_in + 2 * np.pi
                t = np.linspace(t_in, des_t, 100)
                x1, y1 = np.cos(t) * des_r, np.sin(t) * r_in
                plt.plot(x1, y1, color=co)
            else:
                x1 = ori_r * np.cos(ori_t), des_r * np.cos(des_t)
                y1 = ori_r * np.sin(ori_t), des_r * np.sin(des_t)
                plt.plot(x1, y1, color=co)
# ######################### End of plot the map ###############

# ######################### Plot analysis #####################
#  radius of stations: rr
#  angle of stations: tt
# ########## Simplicity
# Number of bends for each line
LineName = ['Disneyland Resort Line', 'East Rail Line', 'Tung Chung Line', 'Island Line',
            'Kwun Tong Line', 'Airport Express', 'Tseung Kwan O Line', 'South Island Line',
            'Tuen Ma Line Phase 1', 'West Rail Line', 'Tsuen Wan Line']
def NumberofBends(MP, Line, line_num, ring_lines, radial_lines):
    NB = []
    for i in range(line_num):
        line_in = np.where(Line == i+1)
        line_index = line_in[0]
        LI = MP[line_index, :]
        LI_in = np.where(LI[:, 1:(ring_lines.size+radial_lines.size + 1)] == 1)
        LI_L = LI_in[1]
        n = 0
        for j in range(LI_L.size):
            if j == 0:
                rl_ID = LI_L[j]
            else:
                if rl_ID != LI_L[j]:
                    if (rl_ID < (ring_lines.size + 1) and rl_ID < (ring_lines.size + 1)) or (
                            rl_ID > (ring_lines.size + 1) and rl_ID > (ring_lines.size + 1)):
                        if j == LI_L.size-1:
                            n = n + 1
                        else:
                            n = n + 2
                    else:
                        n = n + 1
                    rl_ID = LI_L[j]
        NB.append(n)
    return NB
plt.figure()
NB = NumberofBends(MP3, Line, line_num, ring_lines, radial_lines)
plt.bar(range(line_num), NB)
plt.xticks(range(line_num), LineName, rotation=40, fontsize='11')
plt.ylabel('Number of bends')
print('TB:', sum(NB))
print('AB:', sum(NB)/line_num)

# ########## Coherence
# circular used percentage
radial_lines_used = np.zeros((radial_lines.size, 1000))
ring_lines_used = np.zeros((ring_lines.size, 1000))
for i in range(line_num):
    LS = LinkSA[i]
    for j in range(LS.__len__()):  # for each transfer section
        LSS = LS[j]
        for k in range(LSS.__len__()):
            if LSS[k][0] == 1:
                if LSS[k][1] > 0:
                    st = LSS[k][2]
                    et = LSS[k][3]
                    if LSS[k][2] > np.pi:
                        st = LSS[k][2] - 2 * np.pi
                    if LSS[k][2] < -np.pi:
                        st = LSS[k][2] + 2 * np.pi
                    if LSS[k][3] > np.pi:
                        et = LSS[k][3] - 2 * np.pi
                    if LSS[k][3] < -np.pi:
                        et = LSS[k][3] + 2 * np.pi
                    sn = int(st / np.pi * 500)
                    en = int(et / np.pi * 500)
                    ril = np.where(ring_lines == LSS[k][1])
                    if ril[0].shape[0] != 0:
                        if sn * en > 0:
                            if sn < en:
                                for p in range(sn + 499, en + 499, 1):
                                    ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
                            if en == sn:
                                ring_lines_used[ril[0][0], sn + 499] = ring_lines_used[ril[0][0], sn + 499] + 1
                            if sn > en:
                                for p in range(en + 499, sn + 499, 1):
                                    ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
                        else:
                            if sn < en:
                                if sn < - 250:
                                    for p in range(0, sn + 499, 1):
                                        ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
                                    for p in range(en + 499, 999, 1):
                                        ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
                                else:
                                    for p in range(sn + 499, en + 499, 1):
                                        ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
                            if sn > en:
                                if en < -250:
                                    for p in range(0, en + 499, 1):
                                        ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
                                    for p in range(sn + 499, 999, 1):
                                        ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
                                else:
                                    for p in range(en + 499, sn + 499, 1):
                                        ring_lines_used[ril[0][0], p] = ring_lines_used[ril[0][0], p] + 1
            else:
                ral = np.where(abs(radial_lines - LSS[k][1]) < 0.3)
                sn = int(LSS[k][2] * 100) - 1
                en = int(LSS[k][3] * 100) - 1
                if sn < en:
                    for p in range(sn, en, 1):
                        radial_lines_used[ral[0][0], p] = radial_lines_used[ral[0][0], p] + 1
                if en == sn:
                    radial_lines_used[ral[0][0], p] = radial_lines_used[ral[0][0], p] + 1
                if sn > en:
                    for p in range(en, sn, 1):
                        radial_lines_used[ral[0][0], p] = radial_lines_used[ral[0][0], p] + 1
ring_lines_used[ring_lines_used > 0] = 1
CUP = np.sum(ring_lines_used, axis=1)
CUPN = np.sum(CUP) / 1000 / ring_lines.size
print('CUP:', CUPN)
plt.figure()
plt.bar(range(ring_lines.shape[0]), CUP / 1000)
RingLineName = ["RingLine " + str(i + 1) for i in range(ring_lines.shape[0])]
plt.xticks(range(ring_lines.shape[0]), RingLineName, rotation=40, fontsize='11')
plt.ylabel('Occupied percentage')

# ########## Balance
# plot the joint figure
h = sns.JointGrid(x=rr, y=tt)
h.plot_joint(plt.scatter, color=Stationcolor, edgecolor='white')
h.ax_marg_x.hist(rr, color='b', bins=ring_lines)
h.ax_marg_y.hist(tt, color='b', bins=radial_lines, orientation="horizontal")
h.set_axis_labels('radius of ring lines', 'angels of radial line', fontsize=11)
# Length distribution
lenofline = []
lenofallline = []
for i in range(LenInfo.__len__()):
    lenofoneline = []
    for j in range(LenInfo[i].__len__()):
        for k in range(LenInfo[i][j][2]-1):
            lenofoneline.append(LenInfo[i][j][1])
            lenofallline.append(LenInfo[i][j][1])
    lenofline.append(lenofoneline)
length_std = []
length_mean = []
LenIAll = []
plt.figure()
for i in range(line_num):
    me = np.mean(lenofline[i])
    std_deviation = np.std(lenofline[i])
    length_mean.append(me)
    length_std.append(std_deviation)
    plt.errorbar(i, me, yerr=std_deviation, fmt='o', ecolor='r', color='b', elinewidth=2, capsize=4)
plt.xticks(range(11), LineName, rotation=40, fontsize='11')
plt.ylabel('Length')
mean_net = np.mean(lenofallline)
std_net = np.std(lenofallline)
print('Mean_net:', mean_net)
print('Std_net:', std_net)
# Gini coefficient
lengini = copy.deepcopy(lenofallline)
lengini.sort()
cum_lengini = np.cumsum(lengini)
cum_lengini_shares = cum_lengini / sum(lengini)
x = np.linspace(0, 1, lengini.__len__())
B = np.trapz(cum_lengini_shares, x)
A = 0.5 - B
Gini = A/(A+B)
plt.figure()
plt.plot(x, cum_lengini_shares)
plt.xlabel('cumulative share of stations')
plt.ylabel('cumulative share of link length')
print('Gini:', Gini)
# ring /radial line used percentage
used_rr = []
for i in range(line_num):
    LkS = LinkSA[i]
    ring_u = 0
    radial_u = 0
    for j in range(LkS.__len__()):  # for each section between the transfer stations
        LkSp = LkS[j]
        for k in range(LkSp.__len__()):  # line sections between transfers stations
            if LkSp[k][0] == 1:
                if LkSp[k][3] * LkSp[k][2] > 0:
                    ring_u = ring_u + np.abs((LkSp[k][3] - LkSp[k][2]) * LkSp[k][1])
                else:
                    if np.abs(LkSp[k][3]) > np.pi * 0.5:
                        dt = 2 * np.pi - np.abs(LkSp[k][3]) - np.abs(LkSp[k][2])
                    else:
                        dt = np.abs(LkSp[k][3]) + np.abs(LkSp[k][2])
                    ring_u = ring_u + dt * LkSp[k][1]
            else:
                radial_u = radial_u + np.abs(LkSp[k][3] - LkSp[k][2])
    used_rr.append([ring_u, radial_u])
len_ring = 0
len_radial = 0
for i in range(line_num):
    len_ring = len_ring + used_rr[i][0]
    len_radial = len_radial + used_rr[i][1]
RR = len_ring / len_radial
print('RR:', RR)

# # Topographcity
deltal = []
for i in range(tt.__len__()):
    if theta[i] == 0:
        deltal.append(0)
    else:
        dt = abs(theta[i] - tt[i])
        if dt > np.pi:
            deltal.append(2*np.pi - abs(theta[i]) - abs(tt[i]))
        else:
            deltal.append(dt)
al = sum(deltal)/tt.__len__()
print('Topographic:', al)


writer = pd.ExcelWriter('rtheta.xlsx')
rr_data = pd.DataFrame(rr, index=["n"+str(i) for i in range(1, rr.__len__()+1)], columns=["value"])
rr_data.to_excel(writer, 'r', float_format='%.5f')
tt_data = pd.DataFrame(tt, index=["n"+str(i) for i in range(1, rr.__len__()+1)], columns=["value"])
tt_data.to_excel(writer, 't', float_format='%.5f')

writer.save()
writer.close()