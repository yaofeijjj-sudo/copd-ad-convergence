# Flat organ icons for Fig 1, drawn from coordinates here (no third-party artwork). Unit coordinates; cx, cy = centre, s = half-size in mm.
smooth_closed <- function(px, py, n = 400) { k <- length(px); t <- c(seq_len(k), k + 1); s <- seq(1, k + 1, length.out = n)
  data.table(x = spline(t, c(px, px[1]), xout = s, method = "periodic")$y, y = spline(t, c(py, py[1]), xout = s, method = "periodic")$y) }
smooth_open <- function(px, py, n = 60) { t <- seq_along(px); s <- seq(1, length(px), length.out = n)
  data.table(x = spline(t, px, xout = s, method = "natural")$y, y = spline(t, py, xout = s, method = "natural")$y) }
# scalloped contour: push a smooth outline outwards in rounded bumps, as gyri on the brain surface
scallop <- function(d, k, a, from = 0, to = 1) {
  n <- nrow(d); tx <- c(diff(d$x), d$x[1] - d$x[n]); ty <- c(diff(d$y), d$y[1] - d$y[n]); L <- sqrt(tx^2 + ty^2)
  cum <- cumsum(L) / sum(L); on <- cum >= from & cum <= to
  b <- abs(sin(pi * k * (cum - from) / (to - from)))^0.55 * a * on
  sgn <- sign(sum(d$x * c(d$y[-1], d$y[1]) - c(d$x[-1], d$x[1]) * d$y))   # orientation of the contour
  data.table(x = d$x + sgn * (ty / L) * b, y = d$y - sgn * (tx / L) * b)
}
lungs_icon <- function(cx, cy, s) {
  LUNG <- "#E58F85"; LUNG2 <- "#D9766B"; TUBE <- "#B34F4B"; PALE <- "#F6D5D0"
  m <- function(d) d[, .(x = cx + s * x, y = cy + s * y)]
  # viewer's left: right lung with three lobes; viewer's right: left lung with cardiac notch
  R <- smooth_closed(c(-0.30, -0.52, -0.74, -0.90, -1.00, -0.96, -0.74, -0.48, -0.26, -0.13, -0.13, -0.16, -0.20),
                     c(1.02, 0.84, 0.48, 0.04, -0.46, -0.86, -0.90, -0.76, -0.76, -0.56, -0.12, 0.34, 0.78))
  Ln <- smooth_closed(c(0.28, 0.50, 0.72, 0.88, 0.98, 0.94, 0.72, 0.46, 0.30, 0.34, 0.20, 0.13, 0.15, 0.19),
                      c(1.02, 0.84, 0.48, 0.04, -0.46, -0.88, -0.92, -0.80, -0.86, -0.46, -0.20, 0.06, 0.40, 0.80))
  fiss <- list(list(c(-0.17, -0.42, -0.70, -0.96), c(-0.42, -0.28, -0.20, -0.40)), list(c(-0.30, -0.58, -0.84), c(0.20, 0.16, 0.06)),
               list(c(0.17, 0.42, 0.66, 0.92), c(0.14, -0.08, -0.34, -0.66)))
  br <- list(  # main bronchi thick, lobar and segmental branches thinner
    list(c(0, -0.10, -0.26), c(0.50, 0.34, 0.22), 2.2), list(c(0, 0.10, 0.26), c(0.50, 0.34, 0.22), 2.2),
    list(c(-0.26, -0.40, -0.52), c(0.22, 0.42, 0.62), 1.1), list(c(-0.26, -0.44, -0.62), c(0.22, 0.10, 0.02), 1.1), list(c(-0.26, -0.36, -0.50), c(0.22, -0.14, -0.52), 1.3),
    list(c(-0.46, -0.60, -0.66), c(0.48, 0.52, 0.40), 0.6), list(c(-0.40, -0.60, -0.76), c(-0.30, -0.36, -0.56), 0.6), list(c(-0.44, -0.40, -0.34), c(-0.40, -0.56, -0.64), 0.6),
    list(c(0.26, 0.40, 0.50), c(0.22, 0.42, 0.62), 1.1), list(c(0.26, 0.42, 0.64), c(0.22, 0.12, 0.06), 1.1), list(c(0.26, 0.40, 0.60), c(0.22, -0.10, -0.54), 1.3),
    list(c(0.44, 0.58, 0.64), c(0.48, 0.52, 0.38), 0.6), list(c(0.48, 0.66, 0.80), c(-0.26, -0.30, -0.50), 0.6), list(c(0.52, 0.48, 0.42), c(-0.36, -0.56, -0.66), 0.6))
  rings <- seq(0.66, 1.26, length.out = 7)
  c(list(annotate("polygon", x = m(R)$x, y = m(R)$y, fill = LUNG, colour = NA), annotate("polygon", x = m(Ln)$x, y = m(Ln)$y, fill = LUNG, colour = NA)),
    lapply(fiss, function(l) { p <- m(smooth_open(l[[1]], l[[2]], 30)); annotate("path", x = p$x, y = p$y, colour = LUNG2, linewidth = 0.55, lineend = "round") }),
    lapply(br, function(l) { p <- m(smooth_open(l[[1]], l[[2]], 24)); annotate("path", x = p$x, y = p$y, colour = TUBE, linewidth = l[[3]], lineend = "round") }),
    list(annotate("rect", xmin = cx - 0.085 * s, xmax = cx + 0.085 * s, ymin = cy + 0.42 * s, ymax = cy + 1.34 * s, fill = TUBE, colour = NA)),
    lapply(rings, function(r) annotate("segment", x = cx - 0.085 * s, xend = cx + 0.085 * s, y = cy + r * s, yend = cy + r * s, colour = PALE, linewidth = 0.45)))
}
brain_icon <- function(cx, cy, s) {
  CER <- "#7FB2C8"; DK <- "#397A9B"; CB <- "#4F8FAE"; W <- "#FFFFFF"
  m <- function(d) d[, .(x = cx + s * x, y = cy + s * y)]
  base <- smooth_closed(c(-1.02, -0.96, -0.78, -0.50, -0.16, 0.20, 0.54, 0.82, 0.98, 0.98, 0.78, 0.52, 0.30, 0.10, -0.12, -0.34, -0.50, -0.70, -0.90),
                        c(0.02, 0.36, 0.64, 0.82, 0.90, 0.86, 0.72, 0.48, 0.16, -0.12, -0.26, -0.26, -0.42, -0.60, -0.66, -0.56, -0.34, -0.24, -0.16))
  o <- scallop(base, k = 17, a = 0.075)
  cb0 <- smooth_closed(c(0.36, 0.60, 0.84, 0.98, 0.92, 0.70, 0.44, 0.30), c(-0.42, -0.30, -0.30, -0.50, -0.70, -0.82, -0.78, -0.60))
  cbl <- scallop(cb0, k = 11, a = 0.03)
  bs <- smooth_closed(c(0.08, 0.26, 0.40, 0.40, 0.30, 0.16, 0.12, 0.07), c(-0.44, -0.46, -0.66, -1.04, -1.18, -1.16, -0.86, -0.62))
  sul <- list(list(c(-0.56, -0.28, 0.02, 0.32, 0.62, 0.72), c(-0.18, -0.06, -0.04, -0.10, -0.02, 0.10), 1.15),   # lateral sulcus
              list(c(0.02, 0.10, 0.02, 0.10, 0.04), c(0.84, 0.60, 0.40, 0.20, 0.06), 0.95),                       # central sulcus
              list(c(-0.30, -0.36, -0.24, -0.32), c(0.74, 0.50, 0.32, 0.12), 0.8),
              list(c(-0.68, -0.56, -0.68, -0.60), c(0.56, 0.34, 0.20, 0.04), 0.8),
              list(c(0.36, 0.28, 0.40, 0.34), c(0.74, 0.54, 0.36, 0.14), 0.8),
              list(c(0.70, 0.58, 0.68), c(0.48, 0.32, 0.14), 0.8),
              list(c(-0.36, -0.14, 0.08, 0.24), c(-0.34, -0.26, -0.36, -0.28), 0.75),
              list(c(-0.16, -0.40, -0.62), c(0.62, 0.70, 0.52), 0.6), list(c(0.22, 0.38, 0.60), c(0.72, 0.60, 0.60), 0.6))
  fol <- lapply(seq(-0.74, -0.40, length.out = 4), function(k) smooth_open(c(0.34, 0.60, 0.94), c(k + 0.02, k - 0.03, k + 0.07), 30))
  c(list(annotate("polygon", x = m(bs)$x, y = m(bs)$y, fill = DK, colour = NA), annotate("polygon", x = m(cbl)$x, y = m(cbl)$y, fill = CB, colour = NA)),
    lapply(fol, function(p) { p <- m(p); annotate("path", x = p$x, y = p$y, colour = W, linewidth = 0.4, lineend = "round") }),
    list(annotate("polygon", x = m(o)$x, y = m(o)$y, fill = CER, colour = W, linewidth = 0.5)),
    lapply(sul, function(l) { p <- m(smooth_open(l[[1]], l[[2]], 40)); annotate("path", x = p$x, y = p$y, colour = W, linewidth = l[[3]], lineend = "round") }))
}
