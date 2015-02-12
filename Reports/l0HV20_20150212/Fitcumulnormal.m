function [estimates, model] = Fitcumulnormal(xdata, ydata)
% Call fminsearch with a random starting point.
start_point = [0; 1];
model = @cumulnormal;
estimates = fminsearch(model, start_point);
% cumulnormal accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for normcdf-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = cumulnormal(params)
        mu = params(1);
        sigma = params(2);
        FittedCurve = 50+50*normcdf(xdata, mu, sigma);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end