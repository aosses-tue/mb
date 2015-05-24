function G = gauss2d(range_x, range_y, mu_x, mu_y, sigma_x, sigma_y)

g_x = exp(-0.5 * (range_x-mu_x).^2 / sigma_x^2);
g_y = exp(-0.5 * (range_y-mu_y).^2 / sigma_y^2);

G = g_x' * g_y;
