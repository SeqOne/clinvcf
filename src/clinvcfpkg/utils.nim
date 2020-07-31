import logging

var logger* = newConsoleLogger(fmtStr="[$datetime] - $appname - $levelname : ", useStderr=true)
