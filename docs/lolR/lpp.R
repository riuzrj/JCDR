




# 定义LPP函数
lpp <- function(X, k = 2, t = 1) {
  # 计算数据矩阵的大小
  n <- nrow(X)
  d <- ncol(X)
  
  # 计算数据点之间的距离矩阵
  D <- as.matrix(dist(X))
  
  # 计算数据点之间的相似度矩阵
  S <- exp(-D^2/t)
  
  # 计算相似度矩阵的度矩阵
  W <- diag(rowSums(S))
  
  # 计算拉普拉斯矩阵
  L <- W - S
  
  # 计算特征值和特征向量
  eig <- eigen(L, symmetric = TRUE)
  
  # 选择前k个特征向量
  U <- eig$vectors[,1:k]
  
  # 计算降维后的数据矩阵
  Y <- scale(X) %*% U
  
  # 返回结果
  list(Y = Y, U = U, S = S)
}
