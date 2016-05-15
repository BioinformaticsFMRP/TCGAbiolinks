context("Visualize")

test_that("EAbarplot works", {
    Genelist <- c("FN1","COL1A1")
    ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist)
    TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                            GOBPTab = ansEA$ResBP,
                            GOCCTab = ansEA$ResCC,
                            GOMFTab = ansEA$ResMF,
                            PathTab = ansEA$ResPat,
                            nRGTab = Genelist,
                            nBar = 10,
                            filename="a.pdf")
    expect_true(file.exists("a.pdf"))
})

