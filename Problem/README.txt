Example to define a new problem:

prob= smd7mp(1, 2, 1);

upper level evaluation :

[fu, fuc] = prob.evaluate_u(xu, xl)


lower level evaluation
[fl, flc] = prob.evaluate_l(xu, xl)
