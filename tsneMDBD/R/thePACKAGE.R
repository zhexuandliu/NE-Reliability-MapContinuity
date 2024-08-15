########## Singularity Score ##########

#' gradient_compute
#'
#' computes gradients for t-SNE loss function
#'
#' @param Y t-SNE embedding.
#' @param P Similarity matrix in t-SNE.
#' @return Gradients for t-SNE loss function at embedding Y.
#' @examples
#' gradient_compute(Y, P)
#' @export
#' @import ggplot2
#' @import latex2exp
#' @import Rfast
#' @import viridis
#' @export
gradient_compute = function(Y, P){
  Y = c(t(Y))
  n = dim(P)[1]
  n2 = dim(P)[2]
  n_Y = length(Y)
  if ((n!=n2) | (2*n!=n_Y)){
    stop('Input dimensions do not match!')
  }

  YMatmatrix = matrix(rep(c(matrix(matrix(Y, n, 2, byrow = TRUE), 1, 2*n)), n), 2*n, n, byrow = TRUE)
  YMat = matrix(Y, 2, n)
  YnormSq = matrix(colSums(YMat ** 2), n, 1)
  YDiffDistSq = YnormSq %*% matrix(1,1,n) + matrix(1,n,1) %*% t(YnormSq) - 2 * t(YMat) %*% YMat
  YMatmatrix = matrix(rep(c(matrix(matrix(Y, n, 2, byrow = TRUE), 1, 2*n)), n), 2*n, n, byrow = TRUE)
  YDiffDistSq_double = matrix(rep(c(matrix(YDiffDistSq)), each = 2), 2*n, n)
  P_double = matrix(rep(c(matrix(P)), each = 2), 2*n, n)
  deno = sum(1 / (1 + YDiffDistSq)) - n

  I1 = 4 * P_double * (matrix(rep(Y,n),2*n,n) - YMatmatrix) / (1 + YDiffDistSq_double)
  I2 = (-1) * 4 * (matrix(rep(Y,n),2*n,n) - YMatmatrix) / deno / ((1 + YDiffDistSq_double) ^ 2) * sum(P)

  G = rowSums(I1 + I2)

  return(G)
}

#' hessian_compute
#'
#' computes the Hessian matrix for t-SNE loss function
#'
#' @param Y t-SNE embedding.
#' @param P Similarity matrix in t-SNE.
#' @return Hessian matrix for t-SNE loss function at embedding Y.
#' @examples
#' hessian_compute(Y, P)
#' @export
#' @import ggplot2
#' @import latex2exp
#' @import Rfast
#' @import viridis
#' @export
hessian_compute = function(Y, P){
  Y = c(t(Y))
  n = dim(P)[1]
  n2 = dim(P)[2]
  n_Y = length(Y)
  if ((n!=n2) | (2*n!=n_Y)){
    stop('Input dimensions do not match!')
  }

  blockGramMatmatrix = matrix(nrow=2*n, ncol=2*n)
  blockGramMatmatrix2 = matrix(1, nrow=2*n, ncol=2*n)
  indices1 = seq(from=1, to=2*n, by=2)
  indices2 = seq(from=2, to=2*n, by=2)
  blockGramMatmatrix2[indices1, indices2] = 0
  blockGramMatmatrix2[indices2, indices1] = 0
  I2 = matrix(nrow=2*n, ncol=2*n)
  I4 = matrix(nrow=2*n, ncol=2*n)
  I5 = matrix(nrow=2*n, ncol=2*n)

  YMat = matrix(Y, 2, n)
  YnormSq = matrix(colSums(YMat ** 2), n, 1)
  YDiffDistSq = YnormSq %*% matrix(1,1,n) + matrix(1,n,1) %*% t(YnormSq) - 2 * t(YMat) %*% YMat
  deno = sum(1 / (1 + YDiffDistSq)) - n
  YMatmatrix = matrix(rep(c(matrix(matrix(Y, n, 2, byrow = TRUE), 1, 2*n)), n), 2*n, n, byrow = TRUE)
  YDiffDistSq_double = matrix(rep(c(matrix(YDiffDistSq)), each = 2), 2*n, n)

  I1 = ((-4) * P / (1 + YDiffDistSq)) %x% diag(2)
  I3 = 16 * sum(P) / deno**2 * matrix(c(rowSums((matrix(rep(Y,n),2*n,n) - YMatmatrix) / ((1 + YDiffDistSq_double) ^ 2))), 2*n, 1) %*% matrix(c(rowSums((matrix(rep(Y,n),2*n,n) - YMatmatrix) / ((1 + YDiffDistSq_double) ^ 2))), 1, 2*n)
  I4 = (4 * sum(P) / deno / ((1 + YDiffDistSq)**2)) %x% diag(2)

  for (coordinate_idx1 in c(1:2)){
    for (coordinate_idx2 in c(1:2)){
      u = YMat[coordinate_idx1,]
      v = YMat[coordinate_idx2,]
      indices1 = seq(from=coordinate_idx1, to=2*n, by=2)
      indices2 = seq(from=coordinate_idx2, to=2*n, by=2)
      tmp = matrix(u,n,1) %*% matrix(v,1,n)
      blockGramMatmatrix[indices1, indices2] = matrix(rep(diag(tmp),n),n,n) + matrix(rep(diag(tmp),each=n),n,n) - tmp - t(tmp)
      I2[indices1, indices2] = (8 * P) / ((1 + YDiffDistSq)**2) * blockGramMatmatrix[indices1, indices2]
      I5[indices1, indices2] = (-1) * (16 * sum(P)) / deno / ((1 + YDiffDistSq)**3) * blockGramMatmatrix[indices1, indices2]
    }
  }

  H = I1 + I2 + I3 + I4 + I5
  for (coordinate_idx1 in c(1:2)){
    for (coordinate_idx2 in c(1:2)){
      indices1 = seq(from=coordinate_idx1, to=2*n, by=2)
      indices2 = seq(from=coordinate_idx2, to=2*n, by=2)
      diag(H[indices1,indices2])=0
    }
  }
  tmp = apply(array(H, dim=c(2,n,2,n)), c(1,2,3), sum)
  for (i in c(1:n)){
    H[c(2*i-1, 2*i), c(2*i-1, 2*i)] = (-1) * tmp[,i,]
  }

  return (H)
}

#' singularity_score_compute
#'
#' computes the singularity score to diagnose fracture-inducing discontinuity
#'
#' @param Y t-SNE embedding.
#' @param P Similarity matrix in t-SNE.
#' @return Hessian matrix for t-SNE loss function at embedding Y.
#' @examples
#' singularity_score_compute(Y, P)
#' @export
#' @import ggplot2
#' @import latex2exp
#' @import Rfast
#' @import viridis
#' @export
singularity_score_compute = function(Y, P){
  eigen_score_pointwise = function(Hessian, ind){
    h = Hessian[c(2*ind-1, 2*ind), c(2*ind-1, 2*ind)]
    return(1/min(eigen(h)$values))
  }
  H = hessian_compute(Y, P)
  s_score = sapply(c(1:dim(Y)[1]), function(x) eigen_score_pointwise(Hessian = H, x))
  return(s_score)
}

#' plot_singularity_score
#'
#' plot singularity score
#'
#' @param Y t-SNE embedding.
#' @param s_score Singularity score.
#' @param point_size Point size in the plot.
#' @return A plot of t-SNE embedding with singularity score coloring.
#' @examples
#' plot_singularity_score(Y, s_score, point_size = 1)
#' @export
#' @import ggplot2
#' @import latex2exp
#' @import Rfast
#' @import viridis
#' @export
plot_singularity_score = function(Y, s_score, point_size = 1){
  plot_df = data.frame(x = Y[,1], y = Y[,2], s_score = s_score)
  plot_df$quantile_rank = rank(plot_df$s_score) / length(plot_df$s_score)
  p=ggplot() +
    geom_point(data = plot_df,
               aes(x = x, y = y, color = s_score),
               size = point_size) +
    scale_color_viridis(direction = 1,
                        trans = 'log10',
                        name = "Singularity\nScore") +
    xlab('tSNE1') +
    ylab('tSNE2')
  return(p)
}

#' plot_singularity_score_q
#'
#' plot singularity score
#'
#' @param Y t-SNE embedding.
#' @param s_score Singularity score.
#' @param point_size Point size in the plot.
#' @return A plot of t-SNE embedding with the quantile of singularity score as coloring.
#' @examples
#' plot_singularity_score_q(Y, s_score, point_size = 1)
#' @export
#' @import ggplot2
#' @import latex2exp
#' @import Rfast
#' @import viridis
#' @export
plot_singularity_score_q = function(Y, s_score, point_size = 1){
  plot_df = data.frame(x = Y[,1], y = Y[,2], s_score = s_score)
  plot_df$quantile_rank = rank(plot_df$s_score) / length(plot_df$s_score)
  q = c(0,0.25,0.5,0.75,1)
  p=ggplot() +
    geom_point(data = plot_df,
               aes(x = x, y = y, color = quantile_rank),
               size = point_size) +
    scale_color_viridis(direction = 1,
                        breaks = q,
                        labels = as.numeric(sprintf("%.5f", quantile(plot_df$s_score, q))),
                        name = "Singularity\nScore") +
    xlab('tSNE1') +
    ylab('tSNE2')
  return(p)
}

#' plot_singularity_score_q
#'
#' plot singularity score
#'
#' @param perplexity_candidates Candidates for perplexities.
#' @param s_score_mat A matrix of singularity score for each data point at .
#' @param q Point size in the plot.
#' @return A plot of t-SNE embedding with the quantile of singularity score as coloring.
#' @examples
#' plot_s_score_perplexity(perplexity_candidates, s_score_mat, q = 0.95)
#' @export
#' @import ggplot2
#' @import latex2exp
#' @import Rfast
#' @import viridis
#' @import RtsneWithP
#' @export
plot_s_score_perplexity = function(perplexity_candidates, s_score_mat, q = 0.95){
  mean_s_score_vec = apply(s_score_mat, 2, function(s_score) {
    mean(s_score[s_score > quantile(s_score, q)], na.rm = TRUE)
  })
  p = ggplot(data = data.frame(mean = mean_s_score_vec, perplexity = perplexity_candidates)) +
  geom_point(aes(x = perplexity, y = mean))
  return(p)
}

########## Perturbation Score ##########

#' get_loss_pointwise_tsne
#'
#' get surrogate t-SNE loss for specified point
#'
#' @param yy Embedding of the point that we want to calculate loss with.
#' @param Y Embedding matrix of other points.
#' @param Ydist_sq Distance matrix of embedding Y, can be calculated by as.matrix(dist(Y)**2).
#' @param Mat Simalarity matrix (P matrix) of all points (with the point we want to calculate loss with as the first point).
#' @return The surrogate loss for the point yy.
#' @examples
#' get_loss_pointwise_tsne(yy, Y, Mat)
#' @export
#' @import Rfast
#' @import RtsneWithP
#' @export
get_loss_pointwise_tsne = function(yy, Y, Ydist_sq, Mat){
  yy_Y_dist_sq = rowSums((Y - matrix(yy, nrow(Y), ncol(Y), byrow=TRUE))^2)
  yydist_sq = rbind(matrix(yy_Y_dist_sq, nrow=1), as.matrix(Ydist_sq))
  yydist_sq = cbind(matrix(c(0,yy_Y_dist_sq), ncol=1), yydist_sq)
  QMat = (1+yydist_sq)^(-1)
  QMat = QMat / (sum(QMat) - dim(yy)[1])
  eps = 2^(-52)
  PMat = Mat + eps
  QMat = QMat + eps
  return(sum(PMat * log(PMat / QMat)))
}

# perturbation score
# draw contour plot
# run other examples

#' perturbation_score_compute
#'
#' given the perturbation direction and length for specified point, calculate the perturbation score.
#'
#' @param i Calculate the perturbation score for the i-th point.
#' @param X Original Data.
#' @param Y Embedding of original data.
#' @param Ydist_sq Distance matrix of embedding Y, can be calculated by as.matrix(dist(Y)**2).
#' @param perplexity Perplexity parameter to use, should be the same as when calculating the embedding Y.
#' @param pca_result PCA result of original data X, can be calculated by prcomp(X).
#' @param length Length for perturbation.
#' @return Perturbation score for the i-th point in the data.
#' @examples
#' perturbation_score_compute(i, X, Y, Ydist_sq, perplexity, pca_result, length)
#' @export
#' @import RtsneWithP
#' @export
perturbation_score_compute = function(i, X, Y, Ydist_sq, perplexity, pca_result, length){
  PScoreVec = numeric(6)
  dir_id = 1
  if (dim(X)[2] == 2){
    dir_vec = c(
      pca_result$rotation[, 1],-pca_result$rotation[, 1],
      pca_result$rotation[, 2],-pca_result$rotation[, 2]
    )
  }else{
    dir_vec = c(
      pca_result$rotation[, 1],-pca_result$rotation[, 1],
      pca_result$rotation[, 2],-pca_result$rotation[, 2],
      pca_result$rotation[, 3],-pca_result$rotation[, 3]
    )
  }
  for (direction in dir_vec) {
    X_new = rbind(X[i,] + direction * length, X[-i,])
    P_new = Rtsne(X_new,
                  perplexity = perplexity,
                  theta = 0,
                  max_iter = 1,
                  Y_init = Y,
                  check_duplicates = FALSE)
    P_new = P_new$P
    Y_new_init = Y[i,]
    best_value = Inf
    best_result = NULL
    initializations_ind =
      c(i, order(P_new[1,], decreasing = TRUE)[1:2])
    for (init_ind in initializations_ind) {
      if (init_ind < i) {
        Y_new_init = Y[init_ind - 1,]
      } else{
        Y_new_init = Y[init_ind,]
      }
      OPT = optim(Y_new_init, function(x)
        get_loss_pointwise_tsne(x, Y[-i,], Ydist_sq[-i,][,-i], P_new)) # , control = list(maxit = 1000)
      if (OPT$value < best_value) {
        best_value = OPT$value
        best_result = OPT$par
      }
    }
    Y_new = best_result
    PScoreVec[dir_id] = sum((Y_new_init - Y_new) ^ 2) ^ 0.5
    dir_id = dir_id + 1
  }
  return(max(PScoreVec))
}
