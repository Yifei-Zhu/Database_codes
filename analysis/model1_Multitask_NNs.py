import torch
import torch.nn as nn


class Multitask1(nn.Module):
    def __init__(self,nu_features, nu_tasks):
        super(Multitask1,self).__init__()
        self.layer1=nn.Sequential(
                nn.Linear(in_features=nu_features,out_features=1024,bias=True),
                nn.ReLU())
        self.layer2=nn.Sequential(
                nn.Linear(in_features=1024,out_features=512,bias=True),
                nn.ReLU())
        self.layer3=nn.Sequential(
                nn.Linear(in_features=512,out_features=256,bias=True),
                nn.ReLU())
        self.layer4=nn.Sequential(
                nn.Linear(in_features=256,out_features=nu_tasks,bias=True))

    def forward(self,x):
        # x = torch.from_numpy(x).to(torch.float32)
        fc1=self.layer1(x)
        fc2=self.layer2(fc1)
        fc3=self.layer3(fc2)
        output=self.layer4(fc3)
        return output

class Multitask2(nn.Module):
    def __init__(self,nu_features, nu_tasks):
        super(Multitask2,self).__init__()
        self.layer1=nn.Sequential(
                nn.Linear(in_features=nu_features,out_features=512,bias=True),
                nn.ReLU())
        self.layer2=nn.Sequential(
                nn.Linear(in_features=512,out_features=128,bias=True),
                nn.ReLU())
        self.layer3=nn.Sequential(
                nn.Linear(in_features=128,out_features=nu_tasks,bias=True))

    def forward(self,x):
        # x = torch.from_numpy(x).to(torch.float32)
        fc1=self.layer1(x)
        fc2=self.layer2(fc1)
        output=self.layer3(fc2)
        return output


class Multitask3(nn.Module):
    def __init__(self,nu_features, nu_tasks):
        super(Multitask3,self).__init__()
        self.layer1=nn.Sequential(
                nn.Linear(in_features=nu_features,out_features=1024,bias=True),
                nn.ReLU())
        self.layer2=nn.Sequential(
                nn.Linear(in_features=1024,out_features=512,bias=True),
                nn.ReLU())
        self.layer3=nn.Sequential(
                nn.Linear(in_features=512,out_features=256,bias=True),
                nn.ReLU())
        self.layer4=nn.Sequential(
                nn.Linear(in_features=256,out_features=128,bias=True),
                nn.ReLU())
        self.layer5=nn.Sequential(
                nn.Linear(in_features=128,out_features=nu_tasks,bias=True))

    def forward(self,x):
        # x = torch.from_numpy(x).to(torch.float32)
        fc1=self.layer1(x)
        fc2=self.layer2(fc1)
        fc3=self.layer3(fc2)
        fc4=self.layer4(fc3)
        output=self.layer5(fc4)
        return output


class Multitask4(nn.Module):
    def __init__(self,nu_features, nu_tasks):
        super(Multitask4,self).__init__()
        self.layer1=nn.Sequential(
                nn.Linear(in_features=nu_features,out_features=1024,bias=True),
                nn.Sigmoid())
        self.layer2=nn.Sequential(
                nn.Linear(in_features=1024,out_features=512,bias=True),
                nn.Sigmoid())
        self.layer3=nn.Sequential(
                nn.Linear(in_features=512,out_features=256,bias=True),
                nn.Sigmoid())
        self.layer4=nn.Sequential(
                nn.Linear(in_features=256,out_features=128,bias=True),
                nn.Sigmoid())
        self.layer5=nn.Sequential(
                nn.Linear(in_features=128,out_features=nu_tasks,bias=True))

    def forward(self,x):
        # x = torch.from_numpy(x).to(torch.float32)
        fc1=self.layer1(x)
        fc2=self.layer2(fc1)
        fc3=self.layer3(fc2)
        fc4=self.layer4(fc3)
        output=self.layer5(fc4)
        return output


class Multitask5(nn.Module):
    def __init__(self,nu_features, nu_tasks):
        super(Multitask5,self).__init__()
        self.layer1=nn.Sequential(
                nn.Linear(in_features=nu_features,out_features=1024,bias=True),
                nn.ReLU())
        self.layer2=nn.Sequential(
                nn.Linear(in_features=1024,out_features=512,bias=True),
                nn.Sigmoid())
        self.layer3=nn.Sequential(
                nn.Linear(in_features=512,out_features=256,bias=True),
                nn.Sigmoid())
        self.layer4=nn.Sequential(
                nn.Linear(in_features=256,out_features=nu_tasks,bias=True))

    def forward(self,x):
        # x = torch.from_numpy(x).to(torch.float32)
        fc1=self.layer1(x)
        fc2=self.layer2(fc1)
        fc3=self.layer3(fc2)
        output=self.layer4(fc3)
        return output