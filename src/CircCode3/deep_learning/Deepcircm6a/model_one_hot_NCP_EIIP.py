# _*_ coding : utf-8 _*_
import torch
import torch.nn as nn
import torch.nn.functional as F


class CNN51_RNN(nn.Module):
    def __init__(self , HIDDEN_NUM, LAYER_NUM, RNN_DROPOUT, FC_DROPOUT, CELL):
        super(CNN51_RNN ,self).__init__()
        self.basicconv0a = torch.nn.Sequential(
            nn.Conv2d(in_channels=8, out_channels=64, kernel_size=(1, 12), stride=(1,2), padding=(0,2)),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True)
        )# [B, 32, 1, 25]
        self.basicconv0b = torch.nn.Sequential(
            nn.Conv2d(in_channels=64, out_channels=32, kernel_size=(1, 6), stride=(1,2)),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True)
        )# [B, 32, 1, 11]
        self.rnn = BiLSTM_Attention(32 ,HIDDEN_NUM, LAYER_NUM, RNN_DROPOUT)

        self.fc1 = nn.Linear(HIDDEN_NUM * 2, 10)
        self.fc2 = nn.Linear(10, 2)
        self.dropout = nn.Dropout(FC_DROPOUT)

    def forward(self, x):
        x = x.unsqueeze(3).permute(0, 2, 3, 1)
        x = self.basicconv0a(x)
        x = self.basicconv0b(x)
        x = x.squeeze(2).permute(2, 0, 1)  # x > [sequence_len, batch_size, word_vec]
        x = self.rnn(x) # out > [sequence_len, batch_size, num_directions*hidden_size]
        out = self.fc1(x)
        out = self.dropout(out)
        out = F.relu(out)
        out = self.fc2(out)
        return out



class BiLSTM_Attention(nn.Module):
    def __init__(self ,input_size, HIDDEN_NUM, LAYER_NUM, RNN_DROPOUT):
        super(BiLSTM_Attention, self).__init__()
        self.lstm = nn.LSTM(input_size, hidden_size=HIDDEN_NUM, num_layers=LAYER_NUM, bidirectional=True, dropout=RNN_DROPOUT)

    # lstm_output : [batch_size, n_step, HIDDEN_NUM * num_directions(=2)], F matrix
    def attention_net(self, lstm_output, final_state ):
        HIDDEN_NUM = 128
        hidden = final_state.view(-1, HIDDEN_NUM * 2, 3) # hidden : [batch_size, HIDDEN_NUM * num_directions(=2), 3(=n_layer)]
        hidden = torch.mean(hidden, 2).unsqueeze(2)
        attn_weights = torch.bmm(lstm_output, hidden).squeeze(2) # attn_weights : [batch_size, n_step]
        soft_attn_weights = F.softmax(attn_weights, 1)
        context = torch.bmm(lstm_output.transpose(1, 2), soft_attn_weights.unsqueeze(2)).squeeze(2)
        return context # context : [batch_size, HIDDEN_NU

    def forward(self, x):
        input = x

        output, (final_hidden_state, final_cell_state) = self.lstm(input)
        output = output.permute(1, 0, 2) # output : [batch_size, len_seq, HIDDEN_NUM]
        attn_output = self.attention_net(output, final_hidden_state)
        return attn_output

if __name__ == '__main__':
    wordvec_len = 8
    HIDDEN_NUM = 128
    LAYER_NUM = 3
    FC_DROPOUT = 0.5
    RNN_DROPOUT = 0.5
    CELL = 'LSTM'
    X = torch.ones(64,20,8)
    model = CNN51_RNN(HIDDEN_NUM, LAYER_NUM, FC_DROPOUT, RNN_DROPOUT, CELL)
    x = model(X)