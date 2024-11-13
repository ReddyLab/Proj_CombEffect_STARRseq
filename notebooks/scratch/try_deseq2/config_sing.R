FD_RES = "/mount/work/out/scratch"

fun_suppress = function(expr){
    res = suppressWarnings(suppressMessages(eval(expr)))
}