
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split,KFold
from keras.models import Sequential
from keras.layers import Dense, Conv1D, MaxPooling1D, Flatten
from keras.utils import to_categorical
from sklearn.metrics import roc_auc_score, roc_curve

# 加载数据集并进行预处理
data = pd.read_csv('/data/AllDataToBuildModel.csv',index_col=0)
# data.drop(columns=['Row.names'], inplace=True)
X = data.iloc[:, :-1].values
y = data.iloc[:, -1].values


n_splits = 5
kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
auc_scores = []
accuracy_scores = []
for train_idx, test_idx in kf.split(X):
    # 将数据集分为训练集和测试集
    X_train, y_train = X[train_idx], y[train_idx]
    X_test, y_test = X[test_idx], y[test_idx]
    # 将数据集分为训练集和测试集
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
    X_train = X_train.astype('float32')
    X_test = X_test.astype('float32')
    y_train = to_categorical(y_train)
    y_test = to_categorical(y_test)

    # 构建CNN模型
    model = Sequential()
    model.add(Conv1D(filters=32, kernel_size=3, activation='relu', input_shape=(X_train.shape[1], 1)))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(Dense(units=128, activation='relu'))
    model.add(Dense(units=4, activation='softmax'))

    # 编译模型
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy','AUC'])

    # 训练模型
    model.fit(X_train.reshape(X_train.shape[0], X_train.shape[1], 1), y_train, epochs=100, batch_size=32)

    # 在测试集上评估模型的准确性
    loss, accuracy,auc = model.evaluate(X_test.reshape(X_test.shape[0], X_test.shape[1], 1), y_test)
    auc_scores.append(auc)
    accuracy_scores.append(accuracy)

print("交叉验证AUC值:", auc_scores)
print("平均AUC值:", np.mean(auc_scores))
print("交叉验证Accuracy值:", accuracy_scores)
print("平均Accuracy值:", np.mean(accuracy_scores))


# 构建CNN模型
model = Sequential()
model.add(Conv1D(filters=32, kernel_size=3, activation='relu', input_shape=(X.shape[1], 1)))
model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(units=128, activation='relu'))
model.add(Dense(units=4, activation='softmax'))

# 编译模型
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy','AUC'])

# 训练模型
model.fit(X.reshape(X.shape[0], X.shape[1], 1), y, epochs=100, batch_size=32)

# 计算每个基因的贡献度排名
importance = model.layers[0].get_weights()[0].sum(axis=1)
importance_rank = np.argsort(importance)[::-1]
feature_names = data.columns[:-1]
for rank, idx in enumerate(importance_rank):
    print(f"{rank+1}. Feature '{feature_names[idx]}' has importance score {importance[idx]:.4f}")
print('Gene importance rank:', importance_rank)




