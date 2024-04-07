import os
import numpy as np
import pandas as pd
import h5py
import json
import tempfile
import subprocess as sub

from openbabel import pybel
from rdkit import Chem

from sub_code_get_GS_opt_props import main_xyz as get_GS_opt_xyz
from sub_code_get_GS_opt_props import main_import as get_GS_props
from sub_code_get_ES_opt_props import main_import as get_ES_props

def write_molecule_data_to_hdf5(f, ind, smiles1, inchi1, smiles2,inchi2,  properties_ground, properties_excited):
    """
    将一个分子的数据写入HDF5文件

    :param hdf5_file: HDF5文件对象或路径。
    :param cid: 分子的CID。
    :param smiles: 分子的SMILES码。
    :param properties_ground: 分子基态的性质。
    :param properties_excited: 分子激发态的性质。
    """
    # 为当前分子创建一个组
    molecule_group = f.create_group(str(ind))
    molecule_group.attrs['SMILES_PYBEL'] = smiles1
    molecule_group.attrs['INCHI_PYBEL'] = inchi1
    molecule_group.attrs['SMILES_RDKIT'] = smiles2
    molecule_group.attrs['INCHI_RDKIT'] = inchi2

    # 分别为基态和激发态性质创建数据集
    for state_name, properties in [('ground_state', properties_ground), ('excited_state', properties_excited)]:
        state_group = molecule_group.create_group(state_name)
        for prop_name, prop_value in properties.items():

            if isinstance(prop_value, dict):
                prop_json = json.dumps(prop_value)
                dt = h5py.special_dtype(vlen=str)
                # 使用可变长度字符串数据类型来创建数据集
                dset = state_group.create_dataset(prop_name, (1,), dtype=dt, compression='lzf')
                dset[0] = prop_json
            else:
                state_group.create_dataset(prop_name, data=prop_value, compression='lzf')

def xyz_to_SMILESandINCHI_pybel(input_xyz_file):
    molecule = pybel.readfile("xyz", input_xyz_file).__next__()
    smi = molecule.write("smiles")
    inchi = molecule.write("inchi")
    return smi.split()[0], inchi.split()[0]

def xyz_to_SMILESandINCHI_rdkit(xyz_file):
    with tempfile.NamedTemporaryFile(delete=False) as temp_file_mol:
        mol_file = temp_file_mol.name
    sub.run(f'obabel -ixyz {xyz_file} -omol -O {mol_file}', shell=True)
    mol = Chem.MolFromMolFile(mol_file)
    smiles = Chem.MolToSmiles(mol)
    inchi = Chem.MolToInchi(mol)
    os.remove(mol_file)
    return smiles.strip(), inchi.strip()


def get_props_from_log(input_df,hdf5_file):
    err_list = []
    with h5py.File(hdf5_file, 'a') as f:
        for ind in input_df:
            log_file = f'{gau_path_GS}/{ind}/{ind}.log'
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                xyz_file = temp_file.name
            get_GS_opt_xyz(log_file, xyz_file, work_type='opt+freq')
            try:
                smi1, inchi1 = xyz_to_SMILESandINCHI_pybel(xyz_file)
                smi2, inchi2 = xyz_to_SMILESandINCHI_rdkit(xyz_file)
            except:
                err_list.append(ind)
                smi1, inchi1, smi2, inchi2 = '1','1','1','1'
            os.remove(xyz_file)
            log_file_GS=f'{gau_path_GS}/{ind}/{ind}.log'
            log_file_ES=f'{gau_path_ES}/{ind}/{ind}.log'
            properties_ground = get_GS_props(log_file_GS, work_type1)
            properties_excited = get_ES_props(log_file_ES, work_type2)
            write_molecule_data_to_hdf5(f, ind, smi1, inchi1, smi2, inchi2, properties_ground, properties_excited)
    print(err_list)

if __name__ == '__main__':
    # global label, gau_path, state, work_type
    global label, gau_path_GS,  gau_path_ES, work_type1, work_type2

    #label = 'aromatic'
    # label = 'bond_but_nonaromatic'
    label = ''

    work_type1 = 'opt+freq'
    work_type2 = 'sp'


    #main_path = f'/data/home/zhuyf/dataset_work/database/qm9_database'
    main_path = f'/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation'
    #hdf5_file = f'{main_path}/qm9_props.hdf5'
    hdf5_file = f'{main_path}/gdb11_props.hdf5'

    #input_file = f'{main_path}/{label}.csv'
    input_file = f'{main_path}/new_mbkmeans_final_sample_1w.csv'

    #gau_path_GS = f'/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian'
    #gau_path_ES = f'/data/home/zhuyf/dataset_work/database/qm9_database/Gaussian_gau_excited_wb97xd'

    gau_path_GS = f'/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_/gdb11_Gaussian_'
    gau_path_ES = f'/data/home/zhuyf/dataset_work/database/GDB11/split_gdb10/samples_1w_for_calculation/calculation_excited/gdb11_Gaussian_excited'

    #index ='Index'
    index = 'total_index_from_random_sample'

    error_list_path = f'{gau_path_GS}/error_list_{label}.csv'
    error_list= pd.read_csv(error_list_path)['Index']

    input_df = pd.read_csv(input_file)
    input_df = input_df[~input_df[index].isin(error_list)][index]
    #input_df = input_df[~input_df['Index'].isin(error_list)][['Index', 'Smiles_GDB']]

    #input_df = input_df[~input_df['total_index_from_random_sample'].isin(error_list)]['total_index_from_random_sample']



    input_df = [1379421, 2421674, 2182840, 2124462, 467619, 242839, 477046, 966917, 1352003, 228077, 1479014, 419037, 962385, 965587, 566542, 409451, 1019222, 2610063, 2137088, 2607381, 1035344, 242306, 1099739, 1098820, 2024328, 206382, 206390, 127614, 797142, 1735425, 1352009, 2352540, 77082, 348288, 1476710, 206391, 1063956, 834135, 278111, 1190479, 599240, 2831005, 282020, 178241, 233983, 1323519, 1266942, 573248, 270876, 2040669, 130651, 180845, 748420, 1082924, 1355633, 1988494, 694654, 654713, 524788, 226573, 78805, 207615, 438390, 2494035, 1695818, 1488063, 1986213, 1487636, 282194, 81079, 288552, 282247, 1298061, 1558117, 278014, 1323133, 1665867, 1844974, 688350, 688351, 1484834, 476624, 13853, 479989, 129975, 230486, 129992, 387129, 13945, 207663, 235861, 1118889, 933922, 339391, 1191152, 608762, 339498, 605095, 605105, 581363, 1522797, 2919546, 1356760, 79558, 973986, 411781, 1113268, 79587, 1561694, 931210, 1062821, 694412, 213544, 1319575, 1993204, 439242, 285416, 930685, 129402, 789417, 460508, 2019489, 659428, 707683, 702982, 963009, 962992, 790672, 1323436, 685157, 125565, 833922, 1062600, 1081613, 849574, 1998403, 375182, 931457, 935673, 970457, 130771, 180922, 981930, 235561, 1098818, 1996579, 1136246, 410667, 2456808, 623622, 1721174, 1322697, 1721109, 2342832, 1097534, 1044872, 235597, 974373, 1375180, 610635, 1190476, 1375110, 2508019, 1630924, 1629742, 1630948, 1629766, 1375098, 1646803, 1367787, 2336398, 597949, 1328068, 2278136, 425035, 438362, 578838, 573708, 166344, 1743116, 1665156, 1044497, 1045125, 575480, 2261485, 545568, 2378086, 698750, 355455, 1019258, 775344, 281296, 776421, 834129, 835381, 223827, 227343, 207614, 221438, 75949, 227338, 129258, 223826, 835399, 2517427, 554672, 2372257, 946920, 181490, 476745, 2712921, 1065457, 126263, 231288, 289994, 207342, 888360, 165533, 1483029, 243936, 1117059, 232911, 1484802, 208106, 929657, 942284, 836132, 1394638, 1576110, 1314889, 245968, 180837, 220792, 1286178, 694695, 2281323, 101288, 106384, 238745, 1067483, 993367, 494993, 1050407, 1381681, 798685, 894090, 421679, 798677, 1063403, 75241, 75501, 1996127, 2042338, 223892, 134512, 223566, 420617, 1047566, 568749, 221450, 419376, 419672, 278388, 420342, 901064, 1297224, 229269, 1317638, 931056, 2342212, 2112191, 1110529, 1303859, 457123, 1083986, 835149, 1462984, 12792, 74400, 13275, 382598, 289861, 181182, 275501, 287208, 77455, 287253, 170017, 75622, 4351, 286563, 129453, 75681, 438223, 480043, 425779, 411500, 1476621, 1484117, 438283, 236429, 232542, 442959, 228235, 107255, 485229, 168898, 476733, 441688, 2758645, 169958, 223501, 1617876, 1617982, 2595687, 127934, 464787, 2760820, 1472441, 132227, 486243, 580510, 289349, 244741, 421302, 288543, 1395813, 682175, 1629385, 654301, 280587, 166822, 2175358, 1064143, 576189, 476040, 788667, 1709045, 2377206, 2314446, 1652490, 228424, 666796, 2420614, 2420616, 348299, 2306275, 1397522, 933934, 1077821, 381492, 433165, 381024, 208094, 211339, 181460, 850776, 850416, 465838, 381572, 226519, 388646, 1484879, 849589, 834918, 1561765, 974349, 1076230, 1487634, 1487635, 1488061, 207650, 1613252, 2426008, 468360, 468354, 1559930, 1665251, 1381770, 2034111, 1323534, 348060, 2182305, 1036589, 1099777, 178319, 483954, 242822, 90421, 90418, 277991, 170808, 209833, 166422, 240798, 681674, 230484, 134529, 968090, 835461, 129467, 210908, 210911, 206879, 207671, 129439, 981715, 1522363, 1876368, 426472, 974744, 1048979, 974854, 974341, 973976, 974611, 1470255, 2120779, 2135745, 1100235, 446439, 581762, 286546, 79939, 75336, 276810, 170242, 79943, 171225, 1920419, 1006905, 569760, 2374502, 1317625, 345628, 2342333, 79534, 798659, 1829941, 798316, 461454, 893385, 81071, 411550, 106305, 181605, 934429, 1035753, 467445, 441648, 1025896, 276586, 241907, 279059, 276276, 242995, 13852, 287252, 78310, 13944, 2391585, 974001, 84986, 1083460, 968690, 246231, 13748, 247285, 385767, 581842, 208027, 287236, 465645, 476005, 1101626, 486258, 798379, 1630771, 169794, 1052831, 1106284, 1630047, 2600032, 221525, 129278, 833830, 221461, 1077531, 968536, 206855, 1373122, 704879, 893291, 1074650, 1118587, 1303302, 2917710, 232930, 77339, 479419, 439510, 439528, 439519, 424943, 243914, 276060, 901902, 457221, 283384, 2605059, 1098763, 3004508, 289969, 1118965, 429699, 1099009, 1099695, 1099010, 1118706, 1521069, 946872, 1733981, 1561810, 1593694, 2055581, 365865, 1879547, 576419, 382305, 2423873, 385173, 1559129, 208043, 223823, 834809, 230892, 129876, 950776, 970169, 796376, 382304, 1471995, 968693, 967453, 384631, 968951, 382115, 278872, 1463910, 1721432, 384753, 2055284, 1665907, 2032859, 897330, 230880, 798817, 798815, 207286, 134624, 2137086, 2610997, 1630764, 1737682, 230746, 59702, 2124718, 1844744, 1737872, 1558724, 1613328, 1282416, 798852, 1662426, 421897, 1629632, 221618, 2195799, 933879, 1063893, 834816, 1927233, 411861, 430473, 2007408, 2917874, 1713750, 630531, 1069991, 230781, 278150, 569775, 236756, 1737876, 162232, 423951, 430932, 2080260, 2372855, 2585141, 2586660, 2337884, 2024335, 378368, 1291278, 2986991, 1272936, 1657173, 1591624, 1303219, 1915466, 2351916, 1077786, 437073, 1924538, 1586, 2297127, 2391586, 73511, 161199, 1062234, 1401283, 2336099, 1280092, 2597147, 2809935, 2079013, 1075061, 1075577, 1451255, 1075311, 277420, 1078893, 1709044, 1119395, 1463917, 77346, 1078574, 434311, 280872, 1119394, 127711, 421814, 209841, 798381, 2439014, 2020179, 1685419, 1216300, 1129292, 73881, 1048068, 439047, 577569, 438089, 1072477, 1084435, 568297, 1476673, 1062556, 577553, 412180, 577573, 1047625, 162789, 567958, 281298, 1047597, 1048709, 2308080, 2305728, 397335, 2036465, 478984, 478982, 478038, 1321324, 1317103, 347258, 982934, 2423113, 2423602, 2423400, 465403, 654105, 1316659, 1559864, 1377741, 2129458, 337552, 338535, 695525, 241988, 241990, 935162, 2919335, 2774797, 1075502, 1738476, 555089, 558185, 658663, 1986122, 2298396, 1317281, 2307330, 1273206, 1658873, 219848, 705158, 1024539, 273201, 1101934, 1100577, 1101618, 1050360, 1100585, 1101628, 1041746, 1042388, 1046919, 117768, 420556, 420209, 2364019, 1046603, 1062575, 678215, 1091480, 54245, 2065519, 107855, 424912, 441675, 2377628, 438092, 168640, 441656, 420601, 170008, 438276, 475338, 1508326, 2424098, 623636, 1463034, 1063589, 278543, 2644409, 1550485, 588870, 645875, 1998487, 1288596, 1288943, 1297762, 568216, 2130168, 2580547, 1297193, 2175348, 2585288, 485584, 231193, 1323842, 1298330, 1654678, 1589030, 1588621, 1298840, 1654445, 1659127, 1298537, 348349, 446438, 2377005, 688345, 2003779, 163156, 161829, 226478, 420398, 166712, 209836, 276068, 420399, 849564, 117061, 2336700, 1584226, 1630755, 1378147, 1378273, 1630758, 1110329, 2042540, 1042393, 2612037, 820574, 1116729, 1116737, 409254, 397280, 1323046, 1319195, 2659497, 556655, 63823, 586284, 432791, 975689, 1081695, 969211, 1321589, 683109, 2428370, 381710, 902310, 1065067, 2517371, 339458, 337527, 338530, 402895, 2065479, 1317101, 2586759, 656132, 632848, 132229, 1099848, 836519, 2595417, 2595477, 2595332, 412525, 455997, 2458487, 2494384, 2294582, 2363944, 162462, 2372358, 1615792, 1522367, 1522396, 362718, 1658427, 1658535, 1653673, 1653736, 2301914, 2301851, 2059129, 1354647, 484228, 1716329, 2281041, 438382, 1195368, 2971065, 2772889, 467617, 934751, 666639, 1280003, 1267192, 1376751, 2302925, 545048, 46192, 672312, 1268435, 708049, 741780, 576182, 575888, 226879, 243718, 1070214, 1085596, 425637, 434587, 698635, 698910, 281635, 464945, 934987, 130522, 942122, 130528, 849741, 208277, 1110339, 266474, 1508980, 502217, 2165145, 1716922, 1524958, 1477574, 243566, 231017, 465304, 230733, 1484038, 1928264, 1653742, 569753, 426151, 1078900, 970673, 275516, 170818, 77330, 577578, 790178, 280386, 439415, 3690, 40489, 1070166, 439050, 235543, 1048728, 1034956, 479399, 476742, 107208, 434847, 1477552, 1477677, 434814, 168638, 439506, 1037271, 40491, 434409, 569899, 234545, 442941, 2038614, 674426, 2723978, 176798, 577954, 287930, 171248, 280557, 240450, 474509, 579370, 1738444, 1393679, 1717868, 2019491, 223180, 221523, 411484, 419380, 2610972, 955474, 955512, 1081908, 1465028, 1047606, 946900, 282240, 431383, 423305, 1081457, 219445, 1489356, 1081904, 568583, 574790, 234550, 934477, 79965, 43012, 2424327, 229274, 1118482, 1476609, 1023341, 787352, 1476603, 835151, 516003, 516026, 805852, 1456525, 377314, 962335, 376799, 962327, 2296243, 1457161, 2196062, 385593, 102895, 51503, 107796, 208421, 2422044, 1508610, 3007071, 219713, 59713, 547223, 275021, 414357, 967438, 833888, 230163, 2767909, 1662559, 1743182, 1379568, 1321588, 1558319, 162980, 411533, 892667, 432714, 384054, 419792, 411791, 659373, 1617872, 2281569, 2261465, 243378, 2214202, 1588360, 1037500, 835485, 228728, 1069945, 231281, 366975, 2758764, 2250992, 2055761, 897133, 850763, 797923, 78176, 210875, 230156, 2294624, 424791, 485652, 1069977, 276077, 1316701, 206573, 947446, 171456, 281628, 479442, 886278, 1035025, 1025748, 1084486, 882645, 439080, 440772, 581606, 886283, 570239, 576129, 432912, 1778734, 1742769, 1661594, 2030186, 2494418, 210650, 2515729, 2515913, 1380325, 1662550, 388490, 128624, 776113, 1328613, 1330074, 2307356, 1269281, 1328198, 2308002, 1271754, 278672, 1742414, 223259, 830484, 272321, 379498, 1508438, 1100808, 426255, 385789, 1467603, 1117069, 1652089, 657982, 1613342, 1487659, 1323173, 411872, 411886, 882634, 1296779, 630410, 2517395, 484133, 1101609, 1043093, 840891, 1045545, 1044533, 1046173, 445859, 888642, 899914, 1101936, 544458, 127613, 967396, 967363, 797909, 164488, 1034078, 934280, 582461, 1648894, 346272, 1303833, 770297, 2971662, 343811, 2382923, 902604, 378406, 2422496, 901321, 967149, 431394, 1467177, 1083348, 78483, 2424525, 1083243, 423326, 423334, 2493997, 968360, 901315, 2426443, 1467175, 694688, 1588447, 1322145, 2137688, 1395935, 1322367, 2125809, 631122, 385768, 3002257, 1275426, 480016, 1117056, 1742384, 2515606, 1743114, 2378020, 1742385, 1995973, 220683, 75810, 282432, 222921, 206886, 1222968, 955378, 204020, 790067, 1316049, 380648, 968682, 380655, 381033, 966944, 2135880, 367002, 769260, 1980523, 619677, 2428425, 410830, 2772353, 40494, 161699, 485861, 575059, 476613, 1479116, 475455, 485162, 438275, 167024, 1084355, 107846, 581291, 2062412, 424714, 1072603, 1061990, 425617, 1037248, 281301, 1084225, 577562, 40498, 581269, 1083495, 485871, 1063780, 1666448, 432707, 933852, 235596, 486276, 982015, 900411, 2377697, 698693, 2214535, 797984, 129460, 798549, 968101, 1669954, 1670026, 2308472, 2261754, 1557848, 1394796, 2397765, 1875106, 2367955, 1356626, 221586, 2301971, 937577, 835989, 290036, 286560, 290089, 13599, 180767, 969286, 231103, 80352, 1306631, 2772599, 1050083, 948099, 129468, 2305594, 79170, 222576, 2034048, 1717887, 1635175, 621511, 1652145, 769348, 1428853, 892850, 1076012, 1484827, 1048731, 1451348, 1083030, 1490187, 381154, 2673472, 367205, 588966, 2382944, 2586661, 1643016, 576128, 163070, 1453211, 366414, 1733689, 1110753, 1316106, 1033751, 20563, 892894, 411757, 956575, 461684, 898433, 134548, 2774738, 2377684, 2377686, 947442, 467989, 388217, 424454, 431701, 388221, 1467595, 467464, 1466379, 1465529, 466782, 2422962, 388219, 946860, 967192, 948488, 1085583, 103631, 949519, 947454, 903113, 900924, 467442, 431594, 949521, 2423319, 968367, 948132, 1023588, 240099, 612512, 1629839, 1375140, 1330080, 1329975, 1270048, 1109801, 1109803, 931050, 79680, 931013, 1508838, 591763, 2759742, 2759750, 2377165, 2342938, 364486, 163069, 967948, 896174, 612822, 339468, 266388, 963015, 2382776, 1315481, 1618068, 660443, 2342917, 1316238, 231076, 969811, 833819, 222992, 207632, 581481, 2003751, 66353, 2456131, 1667936, 1006571, 2758711, 1652408, 1652396, 930873, 410901, 207648, 127631, 1029675, 246142, 1099507, 838649, 798354, 1100003, 2031387, 366316, 2424440, 2307947, 2306907, 1273152, 1272898, 576180, 221334, 2015865, 1099367, 948103, 1254237, 850760, 381794, 850766, 181455, 211217, 382290, 211337, 209837, 345495, 1471335, 411779, 2299745, 1265064, 1307115, 1647881, 163033, 205010, 968934, 466766, 1049115, 581603, 2421174, 2397305, 75653, 433320, 902143, 132978, 968258, 220500, 223608, 290012, 835790, 244040, 465635, 134623, 244702, 168870, 899156, 244712, 382289, 445820, 235893, 211117, 836002, 798687, 384643, 1896197, 73898, 2421074, 931282, 1036738, 896324, 175074, 455995, 577856, 419102, 221913, 419096, 967439, 223359, 419097, 74551, 420031, 545576, 619676, 366129, 2423487, 203928, 962374, 790361, 17088, 900987, 1844842, 470939, 2428044, 1077938, 2374505, 2374509, 1317624, 1665318, 424830, 239955, 974353, 883593, 226767, 878761, 796271, 275898, 1484878, 59314, 933849, 929378, 1033635, 617137, 1915342, 1550426, 63825, 1972068, 366330, 1356640, 1356654, 1881492, 970164, 421683, 2023758, 2215474, 1571070, 1571065, 2823206, 2597754, 2874356, 928082, 926709, 2336070, 416275, 2398918, 2397700, 2385807, 2277145, 278415, 1656787, 240873, 2515874, 2515665, 1400680, 1743146, 2065518, 108133, 397333, 1267112, 75720, 414417, 1478877, 420565, 420567, 2015823, 1735215, 1083436, 420629, 569161, 224062, 1298346, 346996, 460348, 890932, 2376291, 1464370, 366073, 2195792, 1355659, 258809, 2766650, 1321711, 347757, 426125, 467439, 948481, 385104, 1117045, 468313, 388269, 385480, 472456, 2516342, 1115944, 424938, 969067, 902605, 466789, 983713, 2377367, 679871, 2314368, 599347, 2276426, 1629908, 2505591, 2352649, 364614, 1507691, 281410, 165536, 444774, 465214, 2195949, 1876384, 2031993, 1315126, 1290702, 1316027, 967443, 892986, 221735, 419541, 169851, 892985, 223366, 2130113, 894137, 433335, 1657853, 2064988, 282695, 1049112, 706528, 25869, 221475, 833861, 232512, 127498, 206833, 233991, 835401, 1048326, 1048705, 1470256, 981929, 1082456, 574502, 797147, 1470699, 1076005, 1920263, 2922479, 1778737, 1048706, 850382, 266390, 2280733, 2972216, 1082618, 367054, 1036572, 424452, 983484, 1049430, 568761, 234558, 1465523, 103621, 232482, 221072, 576141, 1489358, 419885, 682364, 572379, 163207, 2091419, 1649245, 1648482, 134525, 230147, 968094, 233996, 2138297, 345713, 37804, 966855, 2423377, 970712, 2065416, 378325, 2423882, 2423374, 2516607, 1742523, 896223, 1078380, 163160, 893419, 1427111, 443056, 485212, 347798, 545565, 545621, 2373608, 2373610, 618549, 2600036, 1099422, 2272069, 2301698, 1462872, 1297023, 108130, 2241777, 429435, 497898, 571701, 1716907, 502237, 439538, 424773, 696520, 2296512, 698593, 600308, 1100802, 2136644, 1303266, 2361579, 2356761, 1370566, 1658429, 1588441, 1658547, 1035275, 1323475, 1508446, 439084, 2007406, 750718, 624113, 624407, 685156, 1398202, 686895, 790705, 2371664, 1322768, 577689, 385734, 969077, 1742694, 439239, 232368, 1006533, 1488414, 228721, 679718, 1400745, 2765703, 3007073, 797992, 207672, 230656, 833889, 221639, 127649, 226472, 221640, 1195366, 222918, 243820, 230745, 798841, 798934, 1377015, 1928883, 382092, 1116394, 1044432, 575529, 850378, 221266, 797885, 277946, 169817, 170787, 850373, 797880, 131867, 81382, 232273, 1118567, 1099924, 950201, 429702, 289964, 381026, 967424, 893391, 381027, 414541, 981719, 829789, 964167, 2122185, 2122108, 3002368, 1741792, 1109849, 1467223, 79714, 967193, 223946, 221815, 419492, 439234, 1484605, 1558127, 421094, 2493296, 2416260, 1617902, 2281330, 659351, 905522, 935196, 1078138, 970170, 1464462, 1471921, 1119017, 835845, 1035054, 2428248, 1085071, 975683, 975976, 975909, 1467771, 1110314, 2428025, 170256, 658007, 661019, 624384, 1323747, 2698511, 708066, 685893, 902176, 388328, 75500, 1077447, 60966, 568770, 968920, 1083496, 461260, 461250, 374954, 1485112, 171233, 798978, 170103, 899142, 206826, 897125, 968034, 419387, 933860, 421276, 933842, 206831, 75479, 1394270, 678572, 1317679, 2748800, 929454, 470756, 925420, 1105225, 232371, 228074, 1508671, 445281, 476006, 1478347, 1101306, 1477875, 204097, 2646146, 388766, 2030916, 2377805, 402892, 1375133, 1613331, 1116604, 1630979, 1322407, 1718174, 433372, 2638703, 1317177, 1099700, 2205521, 2205570, 1082128, 2382822, 706287, 2659821, 2866391, 2582035, 2388149, 748514, 424070, 1528393, 1400738, 3004352, 2758839, 493182, 1928842, 751116, 2356722, 1306825, 275661, 1063951, 2812491, 544395, 2772268, 1330520, 1272564, 599391, 2244476, 1998994, 379226, 2336421, 2335405, 688722, 1931060, 1718106, 1303686, 577575, 1047629, 74443, 1063502, 275513, 962460, 568322, 1085207, 1508841, 229268, 1571009, 724657, 1778873, 646346, 2603761, 2124709, 55962, 1487209, 656450, 656472, 1561325, 1558320, 1476685, 568439, 1484051, 275051, 222564, 206393, 2280769, 2454694, 544999, 970275, 1467219, 969824, 381594, 1084359, 79120, 1085070, 276602, 2765661, 881923, 1091444, 243138, 545388]


    get_props_from_log(input_df, hdf5_file)
    print('Done!')
