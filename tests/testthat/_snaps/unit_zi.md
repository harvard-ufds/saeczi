# printed result is as expected

    Code
      result
    Output
      
      Call:
      unit_zi(samp_dat = samp, pop_dat = pop, lin_formula = lin_formula,      domain_level = "COUNTYFIPS", B = 5, mse_est = TRUE, parallel = FALSE)
      
      Linear Model: 
      - Fixed effects: 
       (Intercept)        tcc16         elev 
      17.787586797  1.238840297 -0.004047397 
      
      - Random effects: 
       Groups     Name        Std.Dev.
       COUNTYFIPS (Intercept) 25.494  
       Residual               68.230  
      
      Logistic Model: 
      - Fixed effects: 
       (Intercept)        tcc16         elev 
      -4.381130205  0.104867536  0.001734101 
      
      - Random effects: 
       Groups     Name        Std.Dev.
       COUNTYFIPS (Intercept) 0.87583 
      

