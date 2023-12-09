
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV, KFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import f1_score, roc_auc_score, accuracy_score

# 读取数据集
data = pd.read_csv('/data/AllDataToBuildModel.csv',index_col=0)
# data.drop(columns=['Row.names'], inplace=True)

# 将基因表达量作为特征矩阵X，将样本的分组情况作为目标变量y
X = data.iloc[:, :-1].values
y = data.iloc[:, -1].values

# 划分训练集、验证集和测试集
X_train_val, X_test, y_train_val, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 划分训练集和验证集（采用五折交叉验证）
kfold = 5
kf = KFold(n_splits=kfold, shuffle=True, random_state=42)

# 定义KNN模型
knn = KNeighborsClassifier()

# 定义参数网格
param_grid = {'n_neighbors': [3, 5, 7]}

# 使用GridSearchCV选择最优参数
grid_search = GridSearchCV(knn, param_grid, cv=kf)
grid_search.fit(X_train_val, y_train_val)

# 获取最优参数
best_params = grid_search.best_params_

# 使用最优参数重新训练模型
knn_best = KNeighborsClassifier(n_neighbors=best_params['n_neighbors'])
knn_best.fit(X_train_val, y_train_val)

# 在测试集上进行预测
y_pred = knn_best.predict(X_test)

# 计算F1分数
f1 = f1_score(y_test, y_pred,average='micro')

# 计算AUC
auc = roc_auc_score(y_test, y_pred,multi_class='ovo')

# 计算准确性
accuracy = accuracy_score(y_test, y_pred)

# 计算各个变量在模型中的贡献度排名
feature_importances = pd.Series(knn_best.feature_importances_, index=X.columns).sort_values(ascending=False)

# 打印评价指标和变量贡献度排名
print('F1 Score:', f1)
print('AUC:', auc)
print('Accuracy:', accuracy)
print('Feature Importances:')
print(feature_importances)
