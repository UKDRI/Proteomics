SETLOCAL ENABLEEXTENSIONS
set selfDir=%~dp0
set methodsDir=%selfDir%..\Methods
dir /a-d /b "%methodsDir%\*.js">"%methodsDir%\methodsList" 
 