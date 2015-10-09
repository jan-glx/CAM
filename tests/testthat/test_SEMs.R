test_that("random polynomial additive SEMs make sense", {
    skip("we do this manually")
    for (i in 1:10) {
        trueDAG <- edges2adj(i=c(1, 2, 3, 2),
                             j=c(3, 3, 4, 4))
        
        sem_object <- random_additive_polynomial_SEM(trueDAG, noise_mean = 1, 
                                                     noise_variance= 0.5,
                                                     intercept_variance=1)
        
        X <- simulate_additive_SEM(sem_object)
        
        pairs(X)
        
        sem_object <- rescale_sem_object(sem_object, 3)
        X <- simulate_additive_SEM(sem_object)
        #png(sprintf("test_ouput/polynomialAdditiveSEM_%40i.png", i))
        pairs(X)
        #dev.off()
        
        nulls <- fitAllOrders(X)
        best_fit <- nulls$bestFits[score== max(nulls$bestFits[,score]),cam]
        best_fit <- best_fit[[sample(length(best_fit),1)]]
        best_fits_resid <- residuals(best_fit)
        
        resid_boot <- CAM:::colwise_resample(best_fits_resid)
        Xboot <- CAM:::fastForward.cam(best_fit, resid_boot)
        pairs(Xboot)
        snd <- fitAllOrders(Xboot)
        best_fit2 <- snd$bestFits[score== max(snd$bestFits[,score]),cam]
        best_fit2 <- best_fit2[[sample(length(best_fit2),1)]]
        logLik(best_fit)
        logLik(best_fit2)
    }
})