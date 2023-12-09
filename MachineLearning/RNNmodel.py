
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd

input_size = 161
hidden_size = 108
num_layers = 1
num_classes = 4

# 加载数据
data = pd.read_csv('/data/DataToBuildModel.csv')
print('Data shape:', data.shape)
# 行是样本，列是基因，最后一列是标签
# 随机选20%的数据作为测试集
test_data = data.sample(frac=0.2, random_state=1)
train_data = data.drop(test_data.index)
test_label = test_data.iloc[:, -1]
test_data = test_data.iloc[:, :-1]
train_label = train_data.iloc[:, -1]
train_data = train_data.iloc[:, :-1]
print('Train data shape:', train_data.shape)
print('Train label shape:', train_label.shape)
print('Label type:',type(train_label))

# 定义RNN模型
class RNNModel(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, output_size):
        super(RNNModel, self).__init__()
        self.hidden_size = hidden_size
        self.rnn = nn.RNN(input_size, hidden_size, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        h0 = torch.zeros(1, x.size(0), self.hidden_size).to(x.device)
        out, _ = self.rnn(x, h0)
        out = self.fc(out[:, -1, :])
        return out

# 定义损失函数和优化器
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = RNNModel(input_size, hidden_size, num_layers, num_classes).to(device)
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 训练模型
num_epochs = 300
for epoch in range(num_epochs):
    for i in range(len(train_data)):
        inputs = torch.tensor(train_data.iloc[i, 1:]).unsqueeze(0).unsqueeze(0).float().to(device)
        labels = torch.tensor(train_label.iloc[i]).long().to(device)
        labels = labels.reshape((1,)).repeat(inputs.shape[1],).to(device)
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
    print('Epoch [{}/{}], Loss: {:.4f}'.format(epoch + 1, num_epochs, loss.item()))

# 计算模型的稳健性
correct = 0
total = 0
with torch.no_grad():
    for i in range(len(test_data)):
        inputs = torch.tensor(test_data.iloc[i, 1:]).unsqueeze(0).unsqueeze(0).float().to(device)
        labels = torch.tensor(test_label.iloc[i]).long().to(device)
        labels = labels.reshape((1,)).repeat(inputs.shape[1],).to(device)
        outputs = model(inputs)
        _, predicted = torch.max(outputs.data, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

accuracy = 100 * correct / total
print('Accuracy of the model on the test set: {:.2f}%'.format(accuracy))

# 给出各个基因的贡献度
weights = model.fc.weight.detach().cpu().numpy()
gene_contributions = np.abs(weights).sum(axis=0) / np.abs(weights).sum()
print('Gene contributions:', gene_contributions)
