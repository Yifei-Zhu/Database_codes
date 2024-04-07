import re

import numpy as np
import pandas as pd
import subprocess as sub

import torch
import torch.nn.functional as F
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from model1_Multitask_NNs import Multitask5 as selectedModel
import torch.utils.data as Data
from sklearn.preprocessing import scale

def readDataFromNpy(descriptor_file):
    descriptor_i = np.load(descriptor_file)
    return np.shape(descriptor_i)[0], np.shape(descriptor_i)[1], descriptor_i

# def processMultitask(task_list):
#     for task in task_list():
#         _, n_dim_y, target = readDataFromNpy()
#     pass

def get_best_device():
    if not torch.cuda.is_available():
        return None
    cmd = "nvidia-smi --query-gpu=memory.total,memory.free,utilization.gpu --format=csv,noheader,nounits"
    result = sub.run(cmd, stdout=sub.PIPE, stderr=sub.PIPE, text=True,shell=True)
    gpus_info = result.stdout.strip().split("\n")
    gpus_used_info = {}
    max_free=0
    for idx, gpu_info in enumerate(gpus_info):
        total_mem, free_mem, gpu_used = gpu_info.split(',')
        if int(free_mem)>max_free:
            best_model=idx

    return idx

def chooseDevice():
    best_device=get_best_device()
    device = torch.device(f"cuda:{best_device}" if best_device != None else "cpu")
    print(f'\nWorking on {device}!\n')
    return device

def setEarlyStop():
    # change the accuracy, loss plot names and model name
    loss_plot_name = 'es_loss'
    acc_plot_name = 'es_accuracy'
    model_name = 'es_model'

def main_train(descriptor_file):
    torch.manual_seed(10)#固定每次初始化模型的权重

    device = chooseDevice()
    training_step = 500
    batch_size = 512
    num_epochs = 10000
    #Earlystop
    best_loss = float('inf')  # 初始化最佳损失为正无穷
    early_stop_patience = 10 # 如果连续3个epoch验证损失没有降低，则停止训练
    counter = 0  # 计数器

    target_file = '/data/home/zhuyf/dataset_work/database/analysis/descriptor/e_homo_lumo_130477samples_2features.npy'
    n_samples, n_features, data = readDataFromNpy(descriptor_file)
    _, n_dim_y, target = readDataFromNpy(target_file)
    # target2 =target1[:, 0] - target1[:, 1]
    # target = np.concatenate((target1, target2.reshape(-1, 1)), axis=1)
    # n_dim_y=np.shape(target)[1]
    # scaled_X = scale(target, axis=0)


    # min_max_scaler = MinMaxScaler()
    # min_max_scaler.fit(data)
    # data = min_max_scaler.transform(data)

    x_train, x_val, y_train, y_val = train_test_split(data,target, test_size=0.2, shuffle=False)

    model = selectedModel(nu_features=n_features, nu_tasks=n_dim_y).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0003)  # 传入 net 的所有参数, 学习率
    loss_func = torch.nn.L1Loss()

    x_train=torch.from_numpy(x_train.astype(np.float32))
    y_train=torch.from_numpy(y_train.astype(np.float32))
    x_val=torch.from_numpy(x_val.astype(np.float32))
    y_val=torch.from_numpy(y_val.astype(np.float32))

    train_data=Data.TensorDataset(x_train,y_train)
    train_data_loader=Data.DataLoader(dataset=train_data,batch_size=64,shuffle=True,num_workers=1)
    val_data=Data.TensorDataset(x_val,y_val)
    val_data_loader=Data.DataLoader(dataset=val_data,batch_size=64,shuffle=True,num_workers=1)

    # early_stopping = EarlyStopping(patience=early_stop_patience, score_function=lambda engine: -engine.state.metrics['val_loss'], trainer=trainer)

    for epoch in range(num_epochs):
        train_data_loader = tqdm(train_data_loader, desc=f'Training Epoch {epoch + 1}/{num_epochs}', leave=True)

        for batch_x, batch_y in train_data_loader:

            batch_x=batch_x.to(device)
            batch_y=batch_y.to(device)

            optimizer.zero_grad()
            outputs = model(batch_x)
            loss = loss_func(outputs, batch_y)
            loss.backward()
            optimizer.step()
            train_data_loader.set_postfix({'Loss': loss.item()}, refresh=True)

        model.eval()  # 设置模型为评估模式
        with torch.no_grad():
            val_loss = 0.0
            for val_batch_inputs, val_batch_labels in val_data_loader:
                val_batch_inputs = val_batch_inputs.to(device)
                val_batch_labels = val_batch_labels.to(device)
                val_outputs = model(val_batch_inputs)
                val_loss += loss_func(val_outputs, val_batch_labels).item()

        val_loss /= len(val_data_loader)
        print(f'validation loss:{val_loss}')
        if val_loss < best_loss:
            best_loss = val_loss
            counter = 0
            torch.save(model.state_dict(), f'{code_path}/Multitask.pth')
        else:
            counter += 1
            if counter >= early_stop_patience:
                print(f'Early stopping after {epoch+1} epochs without improvement on validation set.')
                break

        model.train()
    print(f'\nBest validation loss:{best_loss}\n')
    train_data_loader.close()

def predictNewX():

    X=np.load('/data/home/zhuyf/dataset_work/database/analysis/ecfp.npy')
    # X = np.array([int(bit) for bit in str(ecfp)])
    print(len(X))
    device = chooseDevice()

    model = selectedModel(nu_features=2048, nu_tasks=2).to(device)
    model.load_state_dict(torch.load('Multitask.pth'))

    model.eval()
    X = torch.tensor(X,dtype=torch.float32)
    X = X.to(device)
    with torch.no_grad():
        predictions = model(X)

    # 打印预测结果
    print(predictions)


if __name__ == '__main__':
    global code_path, main_path, gau_path, descriptor_path

    code_path = '/data/home/zhuyf/dataset_work/database/analysis'
    descriptor_path = f'{code_path}/descriptor'

    descriptor_file = f'{descriptor_path}/ECFP4_130477samples_2048features_.npy'
    main_train(descriptor_file)
    # predictNewX()
    print('Done')


