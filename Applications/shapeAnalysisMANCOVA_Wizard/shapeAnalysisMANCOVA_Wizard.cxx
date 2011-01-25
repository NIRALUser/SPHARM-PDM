#include "shapeAnalysisMANCOVA_Wizard.h"
//#define kwsys_stl 
#include <algorithm>
#include <itksys/Process.h>
#include <itksys/SystemTools.hxx>

shapeAnalysisMANCOVA_Wizard::shapeAnalysisMANCOVA_Wizard(std::string infile, QWidget *parent , Qt::WFlags f ) 
: QWidget(parent,f)
{
	setupUi(this);

	if (infile !="nofile") //if csvfile as parameter in the command line
	{
		checkBox_load->setCheckState(Qt::Checked);
			//groupBox_csv->setEnabled(true);
			//frame_del_add_3->setEnabled(true);
			pushButton_load->setEnabled(true);
		QString Infile =QString ((const char*)infile.c_str());
		file_name=Infile;
		infile_ok_display();

	}
	connect(this->checkBox_load, SIGNAL(stateChanged(int)),this, SLOT(mutual_exclusion_load()));
	connect(this->checkBox_NEWCSV, SIGNAL(stateChanged(int)),this, SLOT(mutual_exclusion_new()));

	connect(this->pushButton_load, SIGNAL(clicked()),this, SLOT(openSelectionFileDialog()));	
	connect(this->pushButton_path, SIGNAL(clicked()),this, SLOT(deletePath()));
	connect(this->pushButton_Drow, SIGNAL(clicked()),this, SLOT(deleteRow()));
	connect(this->pushButton_Arow, SIGNAL(clicked()),this, SLOT(addRow()));
	connect(this->pushButton_Dcol, SIGNAL(clicked()),this, SLOT(deleteCol()));
	connect(this->pushButton_Acol, SIGNAL(clicked()),this, SLOT(add_header()));

	connect(this->pushButton_para, SIGNAL(clicked()),this, SLOT(saveParameter()));	
	connect(this->pushButton_folder, SIGNAL(clicked()),this, SLOT(openSelectionFolderDialog()));

	connect(this->pushButton_header, SIGNAL(clicked()),this, SLOT(modifyHeader()));
	connect(tableWidget->horizontalHeader(), SIGNAL(sectionClicked( int)),this, SLOT(openColumnTitleDialog(int)));
	connect(this->checkBox_infile, SIGNAL(stateChanged(int)),this, SLOT(mutual_exclusion_infile(int)));
	connect(this->checkBox_inde, SIGNAL(stateChanged(int)),this, SLOT(mutual_exclusion_inde(int)));
	connect(this->checkBox_gp, SIGNAL(stateChanged(int)),this, SLOT(mutual_exclusion_gp(int)));
	connect(this->checkBox_unselect, SIGNAL(stateChanged(int)),this, SLOT(mutual_exclusion_unseclect(int)));
	connect(this->checkBox_selectscale, SIGNAL(stateChanged(int)),this, SLOT (mutual_exclusion_scale(int)));
	connect(this->tableWidget, SIGNAL(cellClicked( int, int)),this, SLOT(selection_Group_IndeVariabl_DataDialog(int,int)));
	connect(this->checkBox_scale, SIGNAL(stateChanged(int)),this, SLOT(ScaleInTheFile(int)));	
	connect(this->pushButton_save, SIGNAL(clicked()),this, SLOT(saveInfile()));	

	connect(this->checkBox_gptest, SIGNAL(stateChanged(int)),this, SLOT(m_exclu_gptest()));
	connect(this->checkBox_intertest, SIGNAL(stateChanged(int)),this, SLOT(m_exclu_intertest()));	
	connect(this->pushButton_apply, SIGNAL(clicked()),this, SLOT(generate()));

	connect(this, SIGNAL(readyToOpenSlicer()), this, SLOT(openPopUp()));
	connect(this->pushButton_quit, SIGNAL(clicked()), qApp, SLOT(quit()));



}

//*******************************************************************************
//        fonctions
//******************************************************************************




int shapeAnalysisMANCOVA_Wizard::NumberComa()  //test if the file is correct : same number of coma on each line => correct
{    
	int nber_coma_ref,nber_coma;
	nber_coma_ref=0;
	int nber_row;
	int problem = -1;
	int space=0;
	int emptyline=0;
	nber_row=0;
	std::ifstream file(file_name , std::ios::in);

	if (file) {
		std::string line;
       		while(getline(file, line)) 
		{	nber_coma=0;
			space=0;
			for( unsigned int i=0;i<line.length();i++)
			{
				if(line.at(i)==',')
				{	nber_coma++;}
				if(line.at(i)==' ')
				{	space++;}
			}
			if(space !=0 && line.size()==space){emptyline=1;}
			if(problem!=0 && line.size()!=space){
			if(nber_coma_ref==0){nber_coma_ref=nber_coma;}
			else{if(nber_coma!=nber_coma_ref){problem = 0;}}//so problem in the file
			}		
			nber_row++;
			
		}

	}
	nbRow[0]=nber_row;
	if(problem==0){return 0;}	
	else{return (nber_coma_ref);}// so no problem in the file(nber_coma_ref+1=number of columns in the file)

}


void shapeAnalysisMANCOVA_Wizard::InitialisationVectorHeader() //size= number of columns, at the begining fill with 0
{	if(headerVector.size()==0)//at the beginning: fill with 0
	{
		for( int i=0; i<spinBox_col->value();i++)
		{headerVector.push_back(0);}
	}
	else if(spinBox_col->value()>(int)headerVector.size())
	{	for(int i=(int)headerVector.size(); i<spinBox_col->value();i++)
		{headerVector.push_back(0);}//if we add column: new int(=0) in the vector
	}
	else if(spinBox_col->value()<(int)headerVector.size())
	{	for(int i=(int)headerVector.size(); i>spinBox_col->value();i--)
		{headerVector.pop_back();}//if we delete a column, we reduce the size of the vector
	}
}

int shapeAnalysisMANCOVA_Wizard::WhereInfileColumn() //know what is the number of the infile column (return -1 if no infile column)
{
	for(int i=0; i<(int)headerVector.size();i++)
	{if (headerVector[i]==1)return(i);}
	return(-1);
}
int shapeAnalysisMANCOVA_Wizard::WhereScaleColumn()
{
	for(int i=0; i<(int)headerVector.size();i++)
	{if (headerVector[i]==2)return(i);}
	return(-1);
}


void shapeAnalysisMANCOVA_Wizard::infile_ok_display()
{
	nbRow.push_back(0);
	int comment =0;
	InitialisationVectorHeader();
	if(file_name.size()!= 0)
	{lineEdit_path->setText(file_name);
	QString csv= ".csv";
	QString txt=".txt";
	if ((file_name.lastIndexOf(csv) == -1) &&(file_name.lastIndexOf(txt) == -1))
	{
		lineEdit_file->setText(QApplication::translate("MainWindow", " the file is not a .csv/.txt !", 0, QApplication::UnicodeUTF8));
		lineEdit_file->setStyleSheet("background: orange");
	}
	else {
		if(NumberComa()==0)
		{
			lineEdit_file->setText(QApplication::translate("MainWindow", " your csv has a problem", 0, QApplication::UnicodeUTF8));
			lineEdit_file->setStyleSheet("background: orange");
		}
		else
		{
			lineEdit_file->setText(QApplication::translate("MainWindow", " CSV : Ok ", 0, QApplication::UnicodeUTF8));
			lineEdit_file->setStyleSheet("background: lightgreen");
			pushButton_path->setEnabled(true);
			pushButton_load->setEnabled(false);



	
			std::ifstream file(file_name, std::ios::in);//display the file in the grid
			if (file) {
			int col,row,nb;
			col=0;row=0;nb=0;
			std::string line,word;
			spinBox_col->setValue((NumberComa()+1));
			spinBox_data->setValue(nbRow[0]-1);
			ajustCol();
			ajustRow();
			while(getline(file, line)) 
      			{	
				
					for(unsigned int i=0;i<line.length();i++)
					{	
						if ((line.at(i)==',')||(i==line.length()-1))
						{
							QTableWidgetItem* item = new QTableWidgetItem;
							table_display.push_back(item);
							QString qs ( word.c_str() ); //convert the word(string) in QString
							item->setData( 0, qs );
							if(row==0)//first the headers are read
							{
								tableWidget->setHorizontalHeaderItem(col,table_display.back());}
							else
							{tableWidget->setItem( row-1, col, table_display.back());}
							word.clear();
							col++;
						}
						else
						{	word=word+line.at(i);//add a caracter to the word 
							if(line.at(i)=='/' && WhereInfileColumn()==-1 )//if the word contains "/"..the column is the infile column
							{headerVector[col]=1;}
							if (i+1==line.length()-1) {word=word+line.at(i+1);}
						}
						
					}
					if(line.at(0)!='#'){row++;}
					else{comment++;}


			col=0;

			}
			if(WhereInfileColumn()!=-1){lineEditInfileOrScaleDisplay(WhereInfileColumn(),1);}
			spinBox_data->setValue(row-1);
              		file.close(); 
			}
			ajustRow();
			ifEmptyLine();
			pushButton_Drow->setEnabled(true);
			pushButton_Arow->setEnabled(true);
			lineEdit_new_header->setEnabled(true);
			pushButton_Acol->setEnabled(true);
			pushButton_header->setEnabled(true);	
			pushButton_save->setEnabled(true);	
			checkBox_gp->setEnabled(true);	
			checkBox_inde->setEnabled(true);
			checkBox_infile->setEnabled(true);	
			checkBox_scale->setEnabled(true);
			label_infile->setEnabled(true);
			label_group->setEnabled(true);
			//label_inde_var->setEnabled(true);
			checkBox_unselect->setEnabled(true);
			label_unselect->setEnabled(true);
			pushButton_Dcol->setEnabled(true);
			modifyHeaderName.push_back(0);
			tableWidget->setEnabled(true);
			pushButton_Drow->setEnabled(true);
			pushButton_Arow->setEnabled(true);
			InitialisationVectorHeader();
			
		}
	}
}

}


const char* shapeAnalysisMANCOVA_Wizard::QStringToChar (QString myString) //convert QString to char
{
	QByteArray ba=myString.toLatin1();
	const char *myChar=ba.data();
	return(myChar);
}

std::string shapeAnalysisMANCOVA_Wizard::intToString(int myint)//convert int to std::string
{
	std::stringstream mystream;
	mystream<< myint;
	std::string mystr;
	mystr= mystream.str();
	return(mystr);
}

void shapeAnalysisMANCOVA_Wizard::lineEditInfileOrScaleDisplay(int c,int type) // number of the infile column in the QlineEdit and  paint the foreground of the header (idem for scale)
{
	char * char_num_column;
	char_num_column = new char [1];
	strcpy (char_num_column , intToString(c+1).c_str());
	if(type==1)//infile column
	{	lineEdit_infile->setText(QApplication::translate("MainWindow",char_num_column, 0, QApplication::UnicodeUTF8));
		label_infile->setStyleSheet(QString::fromUtf8("color: rgb(39, 74, 255);"));
		headerVector[c]=1;
		paintForeGround(39,74,255,c);
	}
	else//scale
	{	lineEdit_type_sacle->setText(QApplication::translate("MainWindow",char_num_column, 0, QApplication::UnicodeUTF8));
		label_scale_2->setStyleSheet(QString::fromUtf8("color: rgb(190, 49, 213);"));
		headerVector[c]=2;
		paintForeGround(190,49,213,c);
	}

}

void shapeAnalysisMANCOVA_Wizard::paintForeGround(int R,int G,int B,int col)  
{	
	char* header_read = strdup(QStringToChar(tableWidget->horizontalHeaderItem(col)->text()) );
	QTableWidgetItem *header=  new QTableWidgetItem;
	table_display.push_back(header);
	QString str(header_read);
	header->setData( 0, str);
	QBrush brush(QColor(R, G,B, 255));
	brush.setStyle(Qt::SolidPattern);
	header->setForeground(brush);
	tableWidget->setHorizontalHeaderItem(col,table_display.back());
}

void shapeAnalysisMANCOVA_Wizard::deleteLineGoup(int col)  //delete the number col from the vector groupColumn
{
	int where=-1;int i=0;
	do{
		if(groupColumn[i]==col+1)
			{where=i;
			groupColumn.erase (groupColumn.begin()+where);}
		i++;
	}while(where==-1);
}

void shapeAnalysisMANCOVA_Wizard::deleteLineIndependent(int col)  //delete the number col from the vector independentColumn
{
	int where=-1;int i=0;
	do{
		if(independentColumn[i]==col+1)
			{where=i;
			independentColumn.erase (independentColumn.begin()+where);}
		i++;
	}while(where==-1);
	
}

void shapeAnalysisMANCOVA_Wizard::deleteScale(int c) //delete the number c from the infile , no scale in the headerVector
{
	lineEdit_type_sacle->clear();
	paintForeGround(0,0,0,c);
	headerVector[c]=0;
}

void shapeAnalysisMANCOVA_Wizard::displayLineEditIndependent(int c) //diplay the list of the independent variables columns (list is sort and there is , between the numbers)
{
	std::string NumberOfIndeVariableColumn;
	for(unsigned int i=0;i<independentColumn.size();i++)
	{
		char * char_num_column;
		char_num_column = new char [1];
		strcpy (char_num_column , intToString( independentColumn.at(i)).c_str());
		if(i!=0){NumberOfIndeVariableColumn=NumberOfIndeVariableColumn+",";}
		NumberOfIndeVariableColumn=NumberOfIndeVariableColumn+char_num_column;
	}
	char * char_num_column;
	char_num_column = new char [NumberOfIndeVariableColumn.size()+1];
	strcpy (char_num_column, NumberOfIndeVariableColumn.c_str());
	lineEdit_inde_selected->setText(QApplication::translate("MainWindow",char_num_column, 0, QApplication::UnicodeUTF8));
}

void shapeAnalysisMANCOVA_Wizard::displayLineEditGroup(int c)  //diplay the list of the group columns
{
	std::string NumberOfGroupColumn;
	for(unsigned int i=0;i<groupColumn.size();i++)
	{
		std::stringstream num_column_gp;
		num_column_gp<< groupColumn.at(i);;
		std::string str_num_column_gp;
		str_num_column_gp= num_column_gp.str();
		if(i!=0){NumberOfGroupColumn=NumberOfGroupColumn+",";}
		NumberOfGroupColumn=NumberOfGroupColumn+str_num_column_gp;
	}
	char * char_num_column_gp;
	char_num_column_gp = new char [NumberOfGroupColumn.size()+1];
	strcpy (char_num_column_gp, NumberOfGroupColumn.c_str());
	lineEdit_inde_selectedgp->setText(QApplication::translate("MainWindow",char_num_column_gp, 0, QApplication::UnicodeUTF8));
}

void shapeAnalysisMANCOVA_Wizard::ajustRow()  //the Qtablewidget as the same number of columns as the csv
{
	if((spinBox_data->value())>(tableWidget->rowCount()))   
	{
		while((spinBox_data->value())!=(tableWidget->rowCount()))
		{tableWidget->insertRow(tableWidget->rowCount());}
	}

	else if((spinBox_data->value())<(tableWidget->rowCount()) )   
	{
		while((spinBox_data->value())!=(tableWidget->rowCount()))
		{	tableWidget->removeRow(tableWidget->rowCount()-1);}
	}
	
}

void shapeAnalysisMANCOVA_Wizard::ajustCol()  //the Qtablewidget as the same number of rows as the csv
{
	if((spinBox_col->value())>(tableWidget->columnCount()))   
	{
		while((spinBox_col->value())!=(tableWidget->columnCount()))
		{tableWidget->insertColumn(tableWidget->columnCount());}
	}

	else if((spinBox_col->value())<(tableWidget->columnCount()) )   
	{
		while((spinBox_col->value())!=(tableWidget->columnCount()))
		{	tableWidget->removeColumn(tableWidget->columnCount()-1);}
	}
}

void shapeAnalysisMANCOVA_Wizard::saveFile(const char* char_file)  //to store the Qtablewidget in a file
{
	std::ofstream file(char_file, std::ios::out | std::ios::trunc);
        if(file) 
        {
		for(int k=0;k<spinBox_col->value();k++)
		{
			const char *header_read=QStringToChar(tableWidget->horizontalHeaderItem(k)->text());
			if(k!=spinBox_col->value()-1)
			{
				file<<header_read<<",";
			}
			else{
				file<<header_read<<" \n";
			}
		}
		for(int i=0;i<spinBox_data->value();i++)
		{ 	for(int j=0;j<spinBox_col->value();j++)
			{
				tableWidget->setCurrentCell(i,j);
				const char *item_read=QStringToChar(tableWidget->currentItem()->text());
				if(j!=spinBox_col->value()-1)
				{
					file<<item_read<<",";
				}
				else{
					file<<item_read<<" \n";
				}
			}
		}
               file.close();  
        }
}

void shapeAnalysisMANCOVA_Wizard::ifEmptyLine() //if the 1st cell of the last row is empty => empty line at the end =>delete it!
{
	tableWidget->setCurrentCell(spinBox_data->value()-1,0);
	char *item_read=strdup(QStringToChar(tableWidget->currentItem()->text()));
	std::string str(item_read);
		if(str.size()==0)
		{spinBox_data->setValue((spinBox_data->value()-1));
		tableWidget->removeRow((spinBox_data->value()));}
}


void shapeAnalysisMANCOVA_Wizard::setComboBoxGroupColumn()  // dynamic ComboBox
{
	comboBox_testCol->clear();
	for(unsigned int i=0;i<headerVector.size();i++)
	{
		if(headerVector[i]==4)
		{
			comboBox_testCol->addItem(tableWidget->horizontalHeaderItem(i)->text());}
	}
}
void shapeAnalysisMANCOVA_Wizard::setComboBoxIndeVariablesColumn()
{
	comboBox_testCol->clear();

	for(unsigned int i=0;i<headerVector.size();i++)
	{
		if(headerVector[i]==3)
		{ comboBox_testCol->addItem(tableWidget->horizontalHeaderItem(i)->text());}
	}
}



//.*******************************************************************************
//          CSV File visu     File
//.*******************************************************************************


void shapeAnalysisMANCOVA_Wizard::mutual_exclusion_load()
{

	if(checkBox_load->isChecked()==1)
	{	
		if(checkBox_NEWCSV->isChecked()!=1)
		{
			//groupBox_csv->setEnabled(true);
			//frame_del_add_3->setEnabled(true);
			pushButton_load->setEnabled(true);
/*TODO
std::string pathMANCOVA_WizardString;
QString pathMANCOVA_Wizard;
QString pathMANCOVA_Wizard2;
pathMANCOVA_WizardString= itksys::SystemTools::FindProgram("shapeAnalysisMANCOVA_Wizard");
pathMANCOVA_Wizard = pathMANCOVA_WizardString.c_str() ;
pathMANCOVA_Wizard.remove(pathMANCOVA_Wizard.size()-27,27);
QString MANOCAVApicture ="MANOCAVApicture.xpm";
QString Slicerpicture ="Slicer.xpm";
pathMANCOVA_Wizard =pathMANCOVA_Wizard2;
pathMANCOVA_Wizard.append(MANOCAVApicture);
pathMANCOVA_Wizard2.append(Slicerpicture);

label_5->setPixmap(QPixmap(pathMANCOVA_Wizard));
label_6->setPixmap(QPixmap(pathMANCOVA_Wizard));
label_7->setPixmap(QPixmap(pathMANCOVA_Wizard2));*/


		}
		else
		{
			int answer = QMessageBox::question(this, "Warning", "Are you sure you want to load a csv file instead of creating a new one?", QMessageBox::No | QMessageBox::Yes);
			
				if (answer == QMessageBox::Yes)
				{
					checkBox_NEWCSV->setCheckState(Qt::Unchecked);
					//frame_del_add_3->setEnabled(true);
					//groupBox_csv->setEnabled(true);
					//groupBox_new->setEnabled(false);
					pushButton_load->setEnabled(true);
				}
				//else
				if (answer == QMessageBox::No)
				{
					checkBox_NEWCSV->setCheckState(Qt::Checked);
					checkBox_load->setCheckState(Qt::Unchecked);}
		}
	}
}

void shapeAnalysisMANCOVA_Wizard::mutual_exclusion_new()
{
	if(checkBox_NEWCSV->isChecked()==1)
	{
			if(checkBox_load->isChecked()==1)
			{
				int answer = QMessageBox::question(this, "Warning", "Are you sure you want to create a csv file instead of using the one loaded?", QMessageBox::No | QMessageBox::Yes);
				if (answer == QMessageBox::Yes)
    				{
					//groupBox_new->setEnabled(true);
					checkBox_load->setCheckState(Qt::Unchecked);
					//groupBox_csv->setEnabled(false);
					//frame_col2_3->setEnabled(true);
					pushButton_folder->setEnabled(true);
					lineEdit_namecsv->setEnabled(true);
					lineEdit_path->clear();
				}
				
				if (answer == QMessageBox::No)
				{
					checkBox_NEWCSV->setCheckState(Qt::Unchecked);
					checkBox_load->setCheckState(Qt::Checked);}
			}
			else
			{	//groupBox_new->setEnabled(true);
				//frame_col2_3->setEnabled(true);
				pushButton_folder->setEnabled(true);
				lineEdit_namecsv->setEnabled(true);
			}
	}
}





//.*******************************************************************************
//          CSV File visu      Load
//.*******************************************************************************

//find csv file
void shapeAnalysisMANCOVA_Wizard::openSelectionFileDialog()
{

	file_name = QFileDialog::getOpenFileName(this, "select your input", QString());	
	infile_ok_display();
	

}

void shapeAnalysisMANCOVA_Wizard::deletePath()
{
	lineEdit_file->setText(QApplication::translate("MainWindow", " CSV unselected !", 0, QApplication::UnicodeUTF8));
	lineEdit_file->setStyleSheet("background: lightgray");
	lineEdit_path->clear();
	pushButton_path->setEnabled(false);
	pushButton_load->setEnabled(true);
	lineEdit_new_header->setEnabled(false);
	pushButton_Acol->setEnabled(false);
	pushButton_Drow->setEnabled(false);
	pushButton_Arow->setEnabled(false);
	pushButton_Dcol->setEnabled(false);
	if(WhereInfileColumn()!=-1){
		paintForeGround(0,0,0,WhereInfileColumn());
		headerVector[WhereInfileColumn()]=0;}
	table_display.clear();//where all the item are stocked


}

void shapeAnalysisMANCOVA_Wizard::deleteCol()
{
	spinBox_col->setValue((spinBox_col->value()-1));
	tableWidget->removeColumn((spinBox_col->value()));
	if(spinBox_col->value()==4) {pushButton_Dcol->setEnabled(false);}	
}
void shapeAnalysisMANCOVA_Wizard:: addRow()
{
	spinBox_data->setValue((spinBox_data->value()+1));
	tableWidget->insertRow((spinBox_data->value()-1));
}
void shapeAnalysisMANCOVA_Wizard::add_header()
{
	spinBox_col->setValue((spinBox_col->value()+1));
	tableWidget->insertColumn((spinBox_col->value()-1));
	pushButton_Acol->setEnabled(true);

	QTableWidgetItem *header=  new QTableWidgetItem;
	table_display.push_back(header);
	header->setData( 0, lineEdit_new_header->text() );
	tableWidget->setHorizontalHeaderItem((spinBox_col->value()-1),table_display.back());
	pushButton_Dcol->setEnabled(true);

	lineEdit_new_header->clear();
}
void shapeAnalysisMANCOVA_Wizard::deleteRow()
{
	spinBox_data->setValue((spinBox_data->value()-1));
	tableWidget->removeRow((spinBox_data->value()));
}


//.*******************************************************************************
//          CSV File visu      new
//.*******************************************************************************


void shapeAnalysisMANCOVA_Wizard::saveParameter()
{
	InitialisationVectorHeader();//size of the vector= number of columns
	ajustCol();
	ajustRow();
	modifyHeaderName.push_back(0);
	pushButton_header->setEnabled(true);	
	tableWidget->setEnabled(true);	
	checkBox_gp->setEnabled(true);	
	checkBox_inde->setEnabled(true);
	checkBox_infile->setEnabled(true);	
	label_infile->setEnabled(true);
	label_group->setEnabled(true);
	//label_inde_var->setEnabled(true);
	pushButton_save->setEnabled(true);	
	checkBox_unselect->setEnabled(true);
	label_unselect->setEnabled(true);
	checkBox_scale->setEnabled(true);

}

void shapeAnalysisMANCOVA_Wizard::openSelectionFolderDialog()
{
	path = QFileDialog::getExistingDirectory(this);
	lineEdit_path_new->setText(path);
	pushButton_save->setEnabled(true);
}


//.*******************************************************************************
//          CSV File visu      QTableWidget work
//.*******************************************************************************

void shapeAnalysisMANCOVA_Wizard::modifyHeader() //to know we want to change a header(pushbutton) :modifyHeaderName[0]=1, otherwise=0  (size of modifyHeaderName=1)
{
	modifyHeaderName[0]=1;
	checkBox_inde->setChecked(false);
	checkBox_gp->setChecked(false);
	checkBox_infile->setChecked(false);
}

void shapeAnalysisMANCOVA_Wizard::openColumnTitleDialog(int c)  
{	bool ok;
	if (modifyHeaderName[0]==1)//if pushbutton was clicked
	{	
		QString header = QInputDialog::getText("header", "name of the column:", QLineEdit::Normal,
		QString::null, &ok, this );
		if ( ok && !header.isEmpty() ) {
			QTableWidgetItem *header_name=  new QTableWidgetItem;
			table_header.push_back(header_name);
			header_name->setData( 0, header);
			tableWidget->setHorizontalHeaderItem(c,table_header.back());
			if(headerVector[c]==3)//if group column
			{paintForeGround(202,128,35,c);}
			if(headerVector[c]==4)//if inde var column
			{paintForeGround(68,166,60,c);}
			if(headerVector[c]==1)//if infilec column
			{paintForeGround(39,74,255,c);}
			if(headerVector[c]==2)//if scale column
			{paintForeGround(190,49,213,c);}
			modifyHeaderName[0]=0;//we have changed the header
			}
	}
}

void shapeAnalysisMANCOVA_Wizard::mutual_exclusion_infile(int state)
{
	if(state == Qt::Checked) {
		if(checkBox_inde->isChecked()) {checkBox_inde->setCheckState(Qt::Unchecked);}
		if(checkBox_gp->isChecked()) {checkBox_gp->setCheckState(Qt::Unchecked);}
		if(checkBox_selectscale->isChecked()) {checkBox_selectscale->setCheckState(Qt::Unchecked);}
		checkBox_unselect->setCheckState(Qt::Unchecked);
	}
}
void shapeAnalysisMANCOVA_Wizard::mutual_exclusion_inde(int state)
{
	if(state == Qt::Checked) {
		if(checkBox_gp->isChecked()) {checkBox_gp->setCheckState(Qt::Unchecked);}
		if(checkBox_infile->isChecked()) {checkBox_infile->setCheckState(Qt::Unchecked);}
		if(checkBox_selectscale->isChecked()) {checkBox_selectscale->setCheckState(Qt::Unchecked);}
		checkBox_unselect->setCheckState(Qt::Unchecked);
	}
}
void shapeAnalysisMANCOVA_Wizard::mutual_exclusion_gp(int state)
{
	if(state == Qt::Checked) {
		if(checkBox_inde->isChecked()) {checkBox_inde->setCheckState(Qt::Unchecked);}
		if(checkBox_infile->isChecked()) {checkBox_infile->setCheckState(Qt::Unchecked);}
		if(checkBox_selectscale->isChecked()) {checkBox_selectscale->setCheckState(Qt::Unchecked);}
		checkBox_unselect->setCheckState(Qt::Unchecked);

	}
}
void shapeAnalysisMANCOVA_Wizard::mutual_exclusion_unseclect(int state)
{
	if(checkBox_unselect->isChecked()) {
		checkBox_inde->setCheckState(Qt::Unchecked);
		checkBox_gp->setCheckState(Qt::Unchecked);
		checkBox_selectscale->setCheckState(Qt::Unchecked);
		checkBox_infile->setCheckState(Qt::Unchecked);
	}
}
void shapeAnalysisMANCOVA_Wizard::mutual_exclusion_scale(int state)
{
	if(state == Qt::Checked) {
		if(checkBox_gp->isChecked()) {checkBox_gp->setCheckState(Qt::Unchecked);}
		if(checkBox_infile->isChecked()) {checkBox_infile->setCheckState(Qt::Unchecked);}
		if(checkBox_selectscale->isChecked()) {checkBox_inde->setCheckState(Qt::Unchecked);}
		checkBox_unselect->setCheckState(Qt::Unchecked);
	}
}

void shapeAnalysisMANCOVA_Wizard::ScaleInTheFile(int state)
{
	if(state == Qt::Checked) {
		checkBox_selectscale->setEnabled(true);
		label_scale_2->setEnabled(true);
	}
	else
	{
		checkBox_selectscale->setEnabled(false);
		label_scale_2->setStyleSheet(QString::fromUtf8("color: rgb(118,116, 113);"));
		label_scale_2->setEnabled(false);
	}
}

//selection of columns and data
void shapeAnalysisMANCOVA_Wizard::selection_Group_IndeVariabl_DataDialog(int r, int c)  
{	

	if(checkBox_inde->isChecked())
	{ 
		if(headerVector[c]!=1)	// if not the infile column
		{
			if (headerVector[c]==4  ){ //if c was choosen as independent variable column 
				deleteLineGoup(c);
				displayLineEditGroup(c);
				paintForeGround(68,166,60,c);
			}
			if (headerVector[c]==2 ){deleteScale(c);}
			if(headerVector[c]!=3 )//if not already choose as an independent column 
			{
				independentColumn.push_back(c+1);
				sort(independentColumn.begin(), independentColumn.end());
				displayLineEditIndependent(c);
				paintForeGround(202,128,35,c);
				//label_inde_var->setStyleSheet(QString::fromUtf8("color: rgb(202, 128, 35);"));
				headerVector[c]=INDE_COL;
			}
		}
	}


	if(checkBox_gp->isChecked())
	{ 	
		if(headerVector[c]!=1)// if not the infile column
		{
			if (headerVector[c]==3 )//if c was choosen as group column
			{  
				deleteLineIndependent(c);
				displayLineEditIndependent(c);
				paintForeGround(202,128,35,c);}
			if (headerVector[c]==2 ){deleteScale(c);}
			if(headerVector[c]!=4)//if not already choose as a group column
			{	groupColumn.push_back(c+1);
				sort (groupColumn.begin(), groupColumn.end());
				displayLineEditGroup(c);
				paintForeGround(68,166,60,c);
				label_group->setStyleSheet(QString::fromUtf8("color: rgb(68, 166, 60);"));
				headerVector[c]=inde_gp;
			}	
		}
	}

	if(checkBox_infile->isChecked())
	{//id c was choosen as group column or independent independent column it becomes the infile column
		if (headerVector[c]==4)
		{	deleteLineGoup(c);
			displayLineEditGroup(c);
			paintForeGround(68,166,60,c);
		}
		if (headerVector[c]==3)
		{	deleteLineIndependent(c);
			displayLineEditIndependent(c);
			paintForeGround(202,128,35,c);
		}
		if(headerVector[c]==2)
		{deleteScale(c);}
		if (WhereInfileColumn()!=-1){
			paintForeGround(0,0,0,WhereInfileColumn());
			headerVector[WhereInfileColumn()]=0;}
		lineEditInfileOrScaleDisplay(c,1);	
	}

	if(checkBox_selectscale->isChecked())
	{
		if(headerVector[c]!=1)// if not the infile column
		{
			if (headerVector[c]==4)
			{	deleteLineGoup(c);
				displayLineEditGroup(c);
				paintForeGround(68,166,60,c);
			}
			if (headerVector[c]==3)
			{	deleteLineIndependent(c);
				displayLineEditIndependent(c);
				paintForeGround(202,128,35,c);
			}
			if (WhereScaleColumn() !=-1){
				paintForeGround(0,0,0,WhereScaleColumn());
				headerVector[WhereScaleColumn()]=0;}
			lineEditInfileOrScaleDisplay(c,2);
		}
	}
	
	if(checkBox_unselect->isChecked())
	{
		if(headerVector[c]==3)
		{
			deleteLineIndependent(c);
			displayLineEditIndependent(c);
			paintForeGround(0,0,0,c);
			headerVector[c]=0;
		}
		if(headerVector[c]==4)
		{
			deleteLineGoup(c);
			displayLineEditGroup(c);
			paintForeGround(0,0,0,c);
			headerVector[c]=0;
		}
		if(headerVector[c]==2)
		{
			deleteScale(c);
			checkBox_scale->setCheckState(Qt::Unchecked);
		}
	}
	
	if(c ==WhereInfileColumn() && (checkBox_inde->isChecked()!=1 && checkBox_gp->isChecked() !=1 && checkBox_infile->isChecked() !=1  && checkBox_selectscale->isChecked()!=1))
	{ 
		QString text = QFileDialog::getOpenFileName(this, "select your data", QString());
		QTableWidgetItem *item= new QTableWidgetItem;
		path_data.push_back(item);
		item->setData( 0, text );
		tableWidget->setItem ( r, c, path_data.back());
	}
}

void shapeAnalysisMANCOVA_Wizard::saveInfile()
{
	if(WhereInfileColumn()==-1) {QMessageBox::critical(this, "Warning", "choose the column where the paths to the datas are");}
	else{	
		if(checkBox_NEWCSV->isChecked())
		{	QString directory;
			name=lineEdit_namecsv->text();
			if(name.endsWith(".csv")==false){QMessageBox::critical(this, "Warning", "the end of your file name must be '.csv'");}
			else{
				directory=path+'/'+name;
				const char* char_file=QStringToChar(directory);
				saveFile(char_file);}
		}
		else
		{saveFile(file_name);
		} 
	}
}

//.*******************************************************************************
//          Test      (2nd Tab)
//.*******************************************************************************

void shapeAnalysisMANCOVA_Wizard::m_exclu_gptest()
{
	if(checkBox_gptest->isChecked()==1)
	{
		setComboBoxGroupColumn();
		{
			checkBox_intertest->setChecked(false);	
			checkBox_gptest->setChecked(true);
			radioButton_testforcorrelation->setEnabled(false);
			radioButton_negativecorrelation->setEnabled(false);
			radioButton_postivecorrelation->setEnabled(false);
		}
	}		
}
void shapeAnalysisMANCOVA_Wizard::m_exclu_intertest()
{  
	if(checkBox_intertest->isChecked()==1)
	{
		setComboBoxIndeVariablesColumn();
		{
			checkBox_gptest->setChecked(false);	
			checkBox_intertest->setChecked(true);
			radioButton_testforcorrelation->setEnabled(true);
			radioButton_negativecorrelation->setEnabled(true);
			radioButton_postivecorrelation->setEnabled(true);
		}
	}		
}


//pushbutton Apply
void shapeAnalysisMANCOVA_Wizard::generate()
//{
	//emit readyToMANCOVA();
//}
	

//void shapeAnalysisMANCOVA_Wizard::begin()
{	QString qs;
	QProcess *process= new QProcess(this);
	QStringList arguments;


	int numGroup=0; int numInde=0;int infile=0;int testCol=0;int scalecol=0;
	int ComboBoxIndex = comboBox_testCol->currentIndex();
	QString NumColumnGroupTypes=" ";
	QString NumColumnInde = " ";
	QString aide;
	for(unsigned int i=0; i<headerVector.size();i++)
	{
		if(headerVector[i]==3){
			if(checkBox_intertest->isChecked()){
				if(ComboBoxIndex==numInde){testCol=i;}   
				}
			numInde++;
			NumColumnInde=NumColumnInde +" "+ aide.setNum(i);}
		if(headerVector[i]==4){
			if(checkBox_gptest->isChecked()){
				if(ComboBoxIndex==numGroup){testCol=i;}   
				}
			numGroup++;
			NumColumnGroupTypes=NumColumnGroupTypes + " "+aide.setNum(i);
			}
		if(headerVector[i]==1){infile=i;}
		if(headerVector[i]==2){scalecol=i;}
	}



	if(checkBox_load->isChecked())
	{
		arguments.append( file_name);
	}
	else{
		arguments.append( path);
		arguments.append( name);
	}

	if(checkBox_intertest->isChecked()){
		arguments.append("--simpleCorrs");
		arguments.append( "--interactionTest");
		if(radioButton_testforcorrelation->isChecked()){arguments.append("--trendCorrelation");}
		if(radioButton_negativecorrelation->isChecked()){arguments.append("--negativeCorrelation");}
		if(radioButton_postivecorrelation->isChecked()){ arguments.append("--positiveCorrelation");}
	}

	if(radioButton_pillai->isChecked()){ arguments.append("--pillai" );}   
	if(radioButton_hotelling->isChecked()){ arguments.append( "--hotelling" );}
	if(radioButton_wilks->isChecked()){arguments.append( "--wilks" ); }
	if(radioButton_roy->isChecked()){ arguments.append( "--roy" ); }
	
	if(checkBox_debug->isChecked()){ arguments.append( "--debug" );}

	if(checkBox_scale->isChecked()){ 
		arguments.append("--scale");
		qs = QString(intToString(scalecol).c_str());
		arguments.append("--scaleColumn "+qs);
		if(checkBox_computescale->isChecked()){
			arguments.append("--computeScaleFactorFromVolumes");}
	}

	//if(checkBox_zscore->isChecked()){ arguments.append( "-writeZScores" );}


	arguments.append( "--significanceLevel "+lineEdit_pvalue->text());
	
	if(checkBox_KWMinput->isChecked()){arguments.append( "--KWMinput" );}

	qs = QString(intToString(testCol).c_str());
	arguments.append("--testColumn "+qs);

	arguments.append("--columnIndependent "+NumColumnInde);

	qs = QString(intToString(numInde).c_str());
	arguments.append("--numIndependent "+qs);

	arguments.append("--columnGroupTypes "+NumColumnGroupTypes);

	qs = QString(intToString(numGroup).c_str());
	arguments.append("--numGroupTypes "+qs);

	arguments.append("--numPerms "+lineEdit_permu->text());


	qs = QString(intToString(infile).c_str());
	arguments.append("--infileColumn "+qs);

	std::cout<<"-----------------command line shapeAnalysisMANCOVA-----------------"<< std::endl;
	std::cout<<"shapeAnalysisMANCOVA "<< (arguments.join(" ")).toStdString() <<std::endl;
	std::cout<<std::endl;

	std::string pathMANCOVAString;
	QString pathMANCOVA ;
	pathMANCOVAString= itksys::SystemTools::FindProgram("shapeAnalysisMANCOVA");

	if(pathMANCOVAString.empty()==true)
		{ 
			QMessageBox::information(this, "shapeAnalysisMANCOVA", "Select the folder where shapeAnalysisMANCOVA* is saved .");
			pathMANCOVA = QFileDialog::getExistingDirectory(this);
			pathMANCOVA=pathMANCOVA+"/shapeAnalysisMANCOVA";
			QApplication::restoreOverrideCursor();
		}
	else{pathMANCOVA = pathMANCOVAString.c_str() ;}

	QApplication::setOverrideCursor( Qt::WaitCursor );
	process->start ( pathMANCOVA, arguments);

	if( process->waitForFinished(-1)==true){  
		emit readyToOpenSlicer();
	}
}




void shapeAnalysisMANCOVA_Wizard::openPopUp()
{

	QApplication::restoreOverrideCursor();
	int answer = QMessageBox::question(this, "Slicer", "Do you want to open Slicer with the Mrml scene already loaded?", QMessageBox::No | QMessageBox::Yes);

	if (answer == QMessageBox::Yes)
	{
		QString pathSlicer;
		QString pathSlicerCopy;
		std::string pathSlicerString;
		pathSlicerString= itksys::SystemTools::FindProgram("Slicer3");

		//if path not found
		if(pathSlicerString.empty()==true)
		{
			QMessageBox::information(this, "Slicer3", "Select the folder where Slicer3* is saved .");
			pathSlicer = QFileDialog::getExistingDirectory(this);
			//pathSlicer=pathSlicer+"/bin/Slicer3-real";
			
		}
		else{
			std::cout<<" "<<std::endl;
			std::cout<<"path to Slicer"<<pathSlicerString<<std::endl;
		
			//if the Slicer found is in /Slicer/bin/
			std::string key ("bin/Slicer3");
			size_t found;
			found=pathSlicerString.rfind(key);
			if (found!=std::string::npos)
			pathSlicerString.replace (found,key.length(),"Slicer3");

			pathSlicer = pathSlicerString.c_str() ;
			
		}

		QString pathMRML;
		if(checkBox_load->isChecked())
		{
			pathMRML= file_name;
		}
		else{
			pathMRML=path;
			pathMRML.append( name);
		}
		pathMRML.remove(pathMRML.size()-4,4);
		QString end ="_MRMLscene.mrml";
		pathMRML.append(end);
		
		QStringList mrml;mrml.append(QStringToChar(pathMRML));
		mrml << pathMRML;
		QProcess::startDetached(pathSlicer,mrml);
		qApp->quit();
	}


}



//.*******************************************************************************
//         destructor
//.*******************************************************************************
/*
shapeAnalysisMANCOVA_Wizard::~shapeAnalysisMANCOVA_Wizard()
{
	delete Qlineedit_inputfile;
	delete lineEdit_type;
	delete lineEdit_column;
	delete lineEdit_indevar;
	delete lineEdit_outputfile;
//	delete apply ;	
}
*/

