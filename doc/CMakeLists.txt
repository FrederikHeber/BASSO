include (GNUInstallDirs)

set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )

find_package (AsciiDoc)
find_package (DBLATEX)

set (BASSO_DOCDIR ${CMAKE_CURRENT_SOURCE_DIR})
set (BASSO_TITLE ${PROJECT_NAME} ${BASSO_PACKAGE_VERSION})
set (BASSO_STYLESHEET ${BASSO_DOCDIR}/stylesheet.css)
set (BASSO_A2X ${ASCIIDOC_EXECUTABLE} -b docbook -aimagesdir=${CMAKE_CURRENT_SOURCE_DIR} -aversion-label=${PROJECT_NAME} -arevnumber=${BASSO_PACKAGE_VERSION})
set (BASSO_X2D ${DBLATEX_COMMAND})
set (BASSO_A2H ${ASCIIDOC_EXECUTABLE} -b html5 -a latexmath -a stylesheeet=${BASSO_STYLESHEET} -aimagesdir=${CMAKE_CURRENT_SOURCE_DIR} -aversion-label=${PROJECT_NAME} -arevnumber=${BASSO_PACKAGE_VERSION})

macro (add_libbasso_pdf NAME)
add_custom_command (
    OUTPUT ${NAME}.pdf
    COMMAND ${BASSO_A2X} -b docbook -arevdate=$(revdate) -aversion-label=@PACKAGE@ -arevnumber=$(version) -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.xml ${BASSO_DOCDIR}/${NAME}.txt
    COMMAND ${BASSO_X2D} --param=doc.lot.show=figure,table --param=latex.output.revhistory=0 --param=doc.layout="coverpage mainmatter" --param=doc.collab.show=0 --param=paper.type=a4paper -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.pdf ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.xml
    DEPENDS ${BASSO_DOCDIR}/${NAME}*.txt
    DEPENDS ${BASSO_DOCDIR}/examples/*.sh
    DEPENDS ${BASSO_DOCDIR}/listings/*.[ch]pp
    DEPENDS ${BASSO_DOCDIR}/images/*.png
    VERBATIM
)

add_custom_command (
    OUTPUT ${NAME}.html
    COMMAND ${BASSO_A2H} -b html5 -a toc2 -a latexmath -a numbered -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.html ${BASSO_DOCDIR}/${NAME}.txt
    DEPENDS ${BASSO_DOCDIR}/${NAME}*.txt
    DEPENDS ${BASSO_DOCDIR}/examples/*.sh
    DEPENDS ${BASSO_DOCDIR}/listings/*.[ch]pp
    DEPENDS ${BASSO_DOCDIR}/images/*.png
    DEPENDS ${BASSO_STYLESHEET}
    VERBATIM
)

set(BASSO_XML ${BASSO_XML} ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.xml)
set(BASSO_PDF ${BASSO_PDF} ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.pdf)
set(BASSO_HTMLS ${BASSO_HTMLS} ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.html)

install (
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.html
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
)
install (
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.pdf
    DESTINATION ${CMAKE_INSTALL_MANDIR}/${NAME}.pdf
)

endmacro (add_libbasso_pdf)

add_libbasso_pdf (basso)

add_custom_target(run ALL
    DEPENDS basso.pdf basso.html)
