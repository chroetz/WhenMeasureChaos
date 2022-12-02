# Rcpp::cppFunction does not work with slurm. But packages with Rcpp work.
Rcpp::cppFunction('List lorenz63(double t, NumericVector u, NumericVector parms) {
  NumericVector du(3);
  List out(1);
  du[0] = parms[0] * (u[1] - u[0]);
  du[1] = parms[1] * u[0] - u[1] - u[0] * u[2];
  du[2] = u[0] * u[1] - parms[2] * u[2];
  out[0] = du;
  return out;
}')
# lorenz63 <- function(t, u, parms) {
#   du <- c(parms[1] * (u[2] - u[1]),
#           parms[2] * u[1] - u[2] - u[1] * u[3],
#           u[1] * u[2] - parms[3] * u[3])
#   return(list(du))
# }

findAttractor <- function(
    fun, parms,
    u00 = c(1, 1, 20), 
    tMax0 = 100,
    tStep0 = 1e-1
) {
  times0 <- seq(0, tMax0, by = tStep0)
  traj0 <- deSolve::ode(u00, times0, lorenz63, parms)
  return(traj0[nrow(traj0), -1])
}

tMax <- 1e4
tStep <- 1e-3
parms <- c(10, 28, 8 / 3)
u0 <- findAttractor(lorenz63, parms)

# create one long trajectory
times <- seq(0, tMax, by = tStep)
traj <- deSolve::ode(u0, times, lorenz63, parms)
u <- traj[,-1]
n <- nrow(traj)
attractor <- 
  list(
    u = u,
    times = times,
    tStep = tStep,
    n = n)
saveRDS(attractor, "attractorLorenz63.RDS")
