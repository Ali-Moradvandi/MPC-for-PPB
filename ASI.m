%% This function is based on recursive least squere (RLS) to find updated identifed system dynamics
function [theta, P] = ASI(phi,Y,theta,P)
    
    K = P*phi/(1 + phi'*P*phi);

    P = P - K*phi'*P;

    theta = theta + K*(Y - phi'*theta);

end
