# _*_ coding : utf-8 _*_
import torch.nn as nn
import torch.nn.functional as F
from ...deep_learning.end_code_score.seq_load_one_hot import *


class CNN51_RNN(nn.Module):
    def __init__(self, FC_DROPOUT):
        super(CNN51_RNN ,self).__init__()
        self.basicconv0a = torch.nn.Sequential(
            nn.Conv2d(in_channels=4, out_channels=128, kernel_size=(1, 8), stride=(1,2), padding=(0,2)),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True)
        )# [B, 32, 1, 25]
        self.basicconv0b = torch.nn.Sequential(
            nn.Conv2d(in_channels=128, out_channels=32, kernel_size=(1, 6), stride=(1,2)),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True)
        )# [B, 32, 1, 11]
        self.basicconv2a = torch.nn.Sequential(
            nn.Conv2d(in_channels=32, out_channels=32, kernel_size=(1, 4), stride=(1,1)),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True)
        )# [B, 64, 1, 9]

        self.flatten = nn.Flatten()

        self.fc1 = nn.Linear(32 * 20, 10)
        self.fc2 = nn.Linear(10, 2)
        self.dropout = nn.Dropout(FC_DROPOUT)

    def forward(self, x):
        # x > [batch_size, sequence_len, word_vec]
        # print(x.shape)
        x = x.unsqueeze(3).permute(0, 2, 3, 1)
        # print(x.shape)
        x = self.basicconv0a(x)
        x = self.basicconv0b(x)
        x = self.basicconv2a(x) # x > [batch_size, channel(input_size), 1, seq_len]
        #print(x.shape)
        x = self.flatten(x)
        out = self.dropout(self.fc1(x))
        out = F.relu(out)
        out = self.fc2(out)
        return out
