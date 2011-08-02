#ifndef MANCOVAWINDOW_H
#define MANCOVAWINDOW_H

#include "ui_shapeAnalysisMANCOVA_Wizard.h"

#include <itksys/Process.h>

#include <QWidget>
#include <QTextStream>
#include <QFile>
#include <QMessageBox>
#include <QTableWidget>
#include <string.h>
#include <QFileDialog>
#include <QPalette>
#include <iostream>
#include <QPainter>
#include <QHeaderView>
#include <QAbstractItemView>
#include <QBrush>
#include <QTableWidgetItem>
#include <fstream>
#include <QFrame>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <QProcess>
#include <QInputDialog>
#define mypath 1;
#define scale 2;
#define INDE_COL 3;
#define inde_gp 4;

class shapeAnalysisMANCOVA_Wizard : public QWidget, private Ui::MainWindow
{
  Q_OBJECT
public:
  shapeAnalysisMANCOVA_Wizard(std::string infile, QWidget *parent = 0, Qt::WFlags f = 0 );

//	~shapeAnalysisMANCOVA_Wizard();
private slots:

// ****************ON THE1st TAB**************************************************************
  void mutual_exclusion_load(); // the checkbox load is ckecked when the checkbox new is unchecked

  void mutual_exclusion_new(); // the checkbox new is ckecked when the checkbox load is unchecked

// ****************************Load a file***********************************************************
  void openSelectionFileDialog(); // to check if the cvs selected is correct and display it in the table

  void deletePath(); // change the load file //TODO

  void deleteRow(); // pushbutton

  void deleteCol(); // pushbutton

  void addRow(); // pushbutton

  void add_header(); // the name will be the one in the QLineEdit and we validate this with the pushbutton ok

// ********************make a new file****************************************************************
  void saveParameter(); // pushbutton Apply

  void openSelectionFolderDialog(); // select the folder where the file will be save (path to the folder save in the
                                    // QString path)

// ******************working on the Qtablewidget**********************************************************
  void modifyHeader(); // to know we want to change a header(pushbutton) :modifyHeaderName[0]=1, otherwise=0  (size of
                       // modifyHeaderName=1)

  void openColumnTitleDialog(int); // if modifyHeaderName[0]=1 and if we clik on a header => QInputDialog is openning
                                   // and we can write the new name of the header

  void mutual_exclusion_infile(int); // just one checkbox can be checked!

  void mutual_exclusion_inde(int);

  void mutual_exclusion_gp(int);

  void mutual_exclusion_unseclect(int);

  void mutual_exclusion_scale(int);

  void ScaleInTheFile(int);

  void selection_Group_IndeVariabl_DataDialog(int c, int r); // if checkbox is checked, we can select
                                                             // infile/group/inde_var columns or unselect one OR if no
                                                             // checkbox is checked and if you click on a box of the
                                                             // infile column => QFileDialog is openning and we can
                                                             // select the data

  void saveInfile(); // save the QtableWidget as a csv file (as the same place if we have loaded one or in the QString
                     // path with the QString name)

// ******************ON THE 2ND TAB*************************************************
  void m_exclu_gptest();  // just one checkbox can be checked

  void m_exclu_intertest();

  void generate(); // execute shapeAnalysisMANCOVA

  void openPopUp();

signals:
  void readyToOpenSlicer();

  void readyToMANCOVA();

private:
  std::vector<int> modifyHeaderName;  // size of the vector=1. If fill with1->we want to change the name of a header
  QString          pathSlicer;
  QString          file_name; // path of the cvs loaded + name of the file
  QString          path;      // to memorise the path to save the new csv file
  QString          name;      // name of the new cvs file

  std::vector<int> headerVector;  // to know the type of the column : if headerVector[i]=0 : i as no type ; if=1 i is
                                  // the infilecolumn ;if=2 i is a scale column; if= 3 i is a independent variable
                                  // column ; if=4 i is a group column

  std::vector<QTableWidgetItem *> table_display; // where we stock the item of the Qtablewidget

  void infile_ok_display();

  int NumberComa(); // test if the file is correct : same number of coma on each line => correct

  void InitialisationVectorHeader(); // size= number of columns, at the begining fill with 0

  int WhereInfileColumn(); // know what is the number of the infile column (return -1 if no infile column)

  int WhereScaleColumn();

  const char * QStringToChar(QString myString); // convert QString to char

  std::string intToString(int myint); // convert int to std::string

  void lineEditInfileOrScaleDisplay(int c, int type); // number of the infile column in the QlineEdit and  paint the
                                                      // foreground of the header (idem for scale)

  void delete_col(int c);

  void paintForeGround(int R, int G, int B, int col);

  void deleteLineGoup(int col); // delete the number col from the vector groupColumn

  void deleteLineIndependent(int col); // delete the number col from the vector independentColumn

  void deleteScale(int c); // delete the number c from the infile , no scale in the headerVector

  void displayLineEditIndependent(int c); // diplay the list of the independent variables columns (list is sort and
                                          // there is , between the numbers)

  void displayLineEditGroup(int c); // diplay the list of the group columns

  void ajustCol(); // the Qtablewidget as the same number of columns as the csv

  void ajustRow(); // the Qtablewidget as the same number of rows as the csv

  void saveFile(const char* char_file); // to store the Qtablewidget in a file

  void setComboBoxGroupColumn(); // dynamic ComboBox

  void setComboBoxIndeVariablesColumn();

  void ifEmptyLine();

  std::vector<QTableWidgetItem *> table_header; // TODO not usefull
  std::vector<QTableWidgetItem *> path_data;    // TODO not usefull

  std::vector<int> independentColumn; // all the numbers of the independent variable columns
  std::vector<int> groupColumn;       // all the numbers of the group columns
  std::vector<int> nbRow;

  // QString pathSlicer;
  // itksysProcess* m_Process;
  // std::vector<const char*> args;

};

#endif
