Set WshSell = CreateObject("WScript.Shell")
WshSell.Run "cmd /c bodeplots.exe", 0, True
Set WshSell = Nothing