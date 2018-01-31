#include "cdmOptionsWidget.h"

OptionsWidget::OptionsWidget(QMainWindow *_mainwindow, Options *_options)
{

  mainwindow = _mainwindow;
  options = _options;
  name = options->getName();
  description = options->getDescription();
  this->setObjectName("Options");

  connect(this, SIGNAL(updateOptionsWindow()),mainwindow,SLOT(updateOptionsWindow()));

  cfgWindow = new QDialog(this,Qt::Dialog);
  cfgWindow->setModal(true);
  QLabel *nameLabel = new QLabel("Name:");
  nameLineEdit = new QLineEdit;
  QLabel *descriptionLabel = new QLabel("Description:");
  descriptionLineEdit = new QLineEdit;
  QPushButton *finish = new QPushButton("Ok");
  connect(finish,SIGNAL(clicked()),this,SLOT(finish()));
  QPushButton *cancel = new QPushButton("Cancel");
  connect(cancel,SIGNAL(clicked()),this,SLOT(cancel()));
  QGridLayout *cfgWindowlayout = new QGridLayout;
  cfgWindowlayout->addWidget(nameLabel, 0, 0);
  cfgWindowlayout->addWidget(nameLineEdit, 0, 1);
  cfgWindowlayout->addWidget(descriptionLabel, 1, 0);
  cfgWindowlayout->addWidget(descriptionLineEdit, 1, 1);
  cfgWindowlayout->addWidget(finish, 2, 2);
  cfgWindowlayout->addWidget(cancel, 2, 3);
  cfgWindow->setLayout(cfgWindowlayout);
  cfgWindow->setWindowTitle("Save new configuration");

  layout = new QFormLayout(this);

  QHBoxLayout *calculationlayout = new QHBoxLayout();
  executeButton = new QPushButton("Start calculation");
  executeButton->setFixedWidth(150);
  connect(executeButton, SIGNAL(clicked()), this, SLOT(execute()));
  saveButton = new QPushButton("Save configuration");
  saveButton->setFixedWidth(150);
  connect(saveButton, SIGNAL(clicked()), this, SLOT(saveName()));
  calculationlayout->addWidget(executeButton);
  calculationlayout->addWidget(saveButton);

  QHBoxLayout *startnreadlayout = new QHBoxLayout();
  nrigLabel = new QLabel("Calculation options");
  nrig      = new QComboBox();
  nrig->addItems(options->nrigList);
  nrig->setCurrentIndex(options->getNrig());
  nrig->setFixedWidth(150);
  nreadLabel = new QLabel("Read local field from file");
  nread = new QCheckBox(this);
  connect(nread, SIGNAL(stateChanged(int)),this,
	SLOT(nreadCheckBoxStateChanged(int)));
  filerereadLabel = new QLabel("File name:");
  filereread = new QLineEdit(options->getFilereread());
  startnreadlayout->addWidget(filerereadLabel);
  startnreadlayout->addWidget(filereread);
  nmatlabLabel = new QLabel("Do not write mat file");
  nmatlab = new QCheckBox(this);
  connect(nmatlab, SIGNAL(stateChanged(int)),this,
	SLOT(nmatlabCheckBoxStateChanged(int)));
  
  wavelengthLabel = new QLabel("Wavelength (nm)");
  wavelength = new QLineEdit(QString::number(options->getWavelength()));
  wavelength->setFixedWidth(120);
  
  P0Label = new QLabel("P0 (W)");
  P0 = new QLineEdit(QString::number(options->getP0()));
  P0->setFixedWidth(120);

  W0Label = new QLabel("W0 (nm)");
  W0 = new QLineEdit(QString::number(options->getW0()));
  W0->setFixedWidth(120);
  
  beamLabel = new QLabel("Beam");
  beam      = new QComboBox();
  beam->addItems(options->beamList);
  beam->setCurrentIndex(beam->findText(options->getBeam()));
  beam->setFixedWidth(120);
  beam->setCurrentIndex(-1);
  connect(beam , SIGNAL(currentIndexChanged(int)),this,
	SLOT(handleBeamSelectionChanged(int))); 
  beamButton = new QPushButton("Props");
  beamButton->setFixedWidth(70);
  connect(beamButton, SIGNAL(clicked()), this, SLOT(configureBeam()));

  wavemultinumberLabel = new QLabel("Number of plane waves");
  wavemultinumber      = new QSpinBox();
  wavemultinumber->setRange(1,MAX_WAVEMULTI_NUMBER);
  wavemultinumber->setFixedWidth(120);

  objectLabel = new QLabel("Object");
  object      = new QComboBox();
  object->addItems(options->objectList);
  object->setFixedWidth(120);
  object->setCurrentIndex(-1);
  connect(object , SIGNAL(currentIndexChanged(int)),this,
	SLOT(handleObjectSelectionChanged(int)));
  objectButton = new QPushButton("Props");
  objectButton->setFixedWidth(70);
  connect(objectButton, SIGNAL(clicked()), this, SLOT(configureObject()));

  objectnumberLabel = new QLabel("Number of objects");
  objectnumber      = new QSpinBox();
  objectnumber->setRange(1,MAX_OBJECT_NUMBER);
  objectnumber->setFixedWidth(120);

  anisotropyLabel = new QLabel("Anisotropy");
  anisotropy      = new QComboBox();
  anisotropy->addItems(options->anisotropyList);
  anisotropy->setCurrentIndex(anisotropy->findText(options->getAnisotropy()));
  anisotropy->setFixedWidth(120);

  epsilonButton = new QPushButton("Epsilon");
  epsilonButton->setFixedWidth(70);
  connect(epsilonButton, SIGNAL(clicked()), this, SLOT(configureEpsilon()));

  discretizationLabel = new QLabel("Discretization");
  discretization      = new QLineEdit(QString::number(options->getDiscretization()));
  discretization->setFixedWidth(120);

  toleranceLabel = new QLabel("Tolerance");
  tolerance      = new QLineEdit(QString::number(options->getTolerance()));
  tolerance->setFixedWidth(120);

  methodeitLabel = new QLabel("Methode");
  methodeit      = new QComboBox();
  methodeit->addItems(options->methodeitList);
  methodeit->setCurrentIndex(methodeit->findText(options->getMethodeit()));
  methodeit->setFixedWidth(120);

  polarizabilityLabel = new QLabel("Polarizability");
  polarizability      = new QComboBox();
  polarizability->addItems(options->polarizabilityList);
  polarizability->setCurrentIndex(polarizability->findText(options->getPolarizability()));
  polarizability->setFixedWidth(120);

  quadLabel = new QLabel("Quad");
  quad      = new QComboBox();
  quad->addItems(options->quadList);
  quad->setCurrentIndex(quad->findText(QString::number(options->getQuad())));
  quad->setFixedWidth(120);


  nfft2dLabel = new QLabel("FFT");
  nfft2d      = new QComboBox();
  nfft2d->addItems(options->nfft2dList);
  nfft2d->setCurrentIndex(nfft2d->findText(QString::number(options->getnfft2d())));
  nfft2d->setFixedWidth(120);
  
  QGridLayout *studyfarfieldlayout = new QGridLayout();
  emptycrosssectionLabel = new QLabel(" ");
  crosssectionLabel = new QLabel("Cross section");
  crosssection = new QCheckBox(this);
  connect(crosssection, SIGNAL(stateChanged(int)),this,
	SLOT(crosssectionCheckBoxStateChanged(int)));
  emptycrosssectionpoyntingLabel = new QLabel(" ");
  crosssectionpoyntingLabel = new QLabel("Cross section + Poynting");
  crosssectionpoynting = new QCheckBox(this);
  connect(crosssectionpoynting, SIGNAL(stateChanged(int)),this,
	SLOT(crosssectionpoyntingCheckBoxStateChanged(int)));
  quickdiffractLabel = new QLabel("Quick computation");
  quickdiffract = new QCheckBox(this);
  ntheta = new QLineEdit(QString::number(options->getNtheta()));
  ntheta->setFixedWidth(40);
  nthetaLabel = new QLabel("Ntheta:");
  nphi = new QLineEdit(QString::number(options->getNphi()));
  nphi->setFixedWidth(40);
  nphiLabel = new QLabel("Nphi:");
  emptynenergieLabel = new QLabel(" ");
  nenergieLabel = new QLabel("Emissivity");
  nenergie = new QCheckBox(this);
  connect(nenergie, SIGNAL(stateChanged(int)),this,
	SLOT(nenergieCheckBoxStateChanged(int)));
  emptymicroscopyLabel = new QLabel(" ");
  microscopyLabel = new QLabel("Microscopy");
  microscopy = new QCheckBox(this);
  connect(microscopy, SIGNAL(stateChanged(int)),this,
	SLOT(microscopyCheckBoxStateChanged(int)));
  microscopyFFTLabel = new QLabel("Quick computation");
  microscopyFFT = new QCheckBox(this);
  na = new QLineEdit(QString::number(options->getNA()));
  na->setFixedWidth(40);
  naLabel = new QLabel("Numerical aperture [0,1]:");
  gross = new QLineEdit(QString::number(options->getGross()));
  gross->setFixedWidth(40);
  grossLabel = new QLabel("Magnification:");
  studyfarfieldlayout->addWidget(emptycrosssectionLabel,0,0,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssectionLabel,0,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssection,0,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(emptycrosssectionpoyntingLabel,1,0,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssectionpoyntingLabel,1,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssectionpoynting,1,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(quickdiffractLabel,1,3,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(quickdiffract,1,4,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nthetaLabel,2,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(ntheta,2,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nphiLabel,3,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nphi,3,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(emptynenergieLabel,4,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nenergieLabel,4,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nenergie,4,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(emptymicroscopyLabel,5,0,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(microscopyLabel,5,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(microscopy,5,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(microscopyFFTLabel,5,3,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(microscopyFFT,5,4,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(naLabel,6,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(na,6,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(grossLabel,7,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(gross,7,2,Qt::AlignLeft);
   
  QGridLayout *studyforcelayout = new QGridLayout();
  emptyopticalforceLabel = new QLabel(" ");
  opticalforceLabel = new QLabel("Optical force");
  opticalforce = new QCheckBox(this);
  connect(opticalforce, SIGNAL(stateChanged(int)),this,
	SLOT(opticalforceCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticalforceLabel,0,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforceLabel,0,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforce,0,2,Qt::AlignLeft);
  emptyopticalforcedensityLabel = new QLabel(" ");
  opticalforcedensityLabel = new QLabel("Optical force density");
  opticalforcedensity = new QCheckBox(this);
  connect(opticalforcedensity, SIGNAL(stateChanged(int)),this,
	SLOT(opticalforcedensityCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticalforcedensityLabel,1,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforcedensityLabel,1,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforcedensity,1,2,Qt::AlignLeft);
  emptyopticaltorqueLabel = new QLabel(" ");
  opticaltorqueLabel = new QLabel("Optical torque");
  opticaltorque = new QCheckBox(this);
  connect(opticaltorque, SIGNAL(stateChanged(int)),this,
	SLOT(opticaltorqueCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticaltorqueLabel,2,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorqueLabel,2,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorque,2,2,Qt::AlignLeft);
  emptyopticaltorquedensityLabel = new QLabel(" ");
  opticaltorquedensityLabel = new QLabel("Optical torque density");
  opticaltorquedensity = new QCheckBox(this);
  connect(opticaltorquedensity, SIGNAL(stateChanged(int)),this,
	SLOT(opticaltorquedensityCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticaltorquedensityLabel,3,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorquedensityLabel,3,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorquedensity,3,2,Qt::AlignLeft);

  QGridLayout *studynearfieldlayout = new QGridLayout();
  emptylocalfieldLabel = new QLabel(" ");
  localfieldLabel = new QLabel("Local field");
  localfield = new QCheckBox();
  connect(localfield, SIGNAL(stateChanged(int)),this,
	SLOT(localfieldCheckBoxStateChanged(int)));
  studynearfieldlayout->addWidget(emptylocalfieldLabel,0,0,Qt::AlignLeft);
  studynearfieldlayout->addWidget(localfieldLabel,0,1,Qt::AlignLeft);
  studynearfieldlayout->addWidget(localfield,0,2,Qt::AlignLeft);
  emptymacroscopicfieldLabel = new QLabel(" ");
  macroscopicfieldLabel = new QLabel("Macroscopic field");
  macroscopicfield = new QCheckBox();
  connect(macroscopicfield, SIGNAL(stateChanged(int)),this,
	SLOT(macroscopicfieldCheckBoxStateChanged(int)));
  studynearfieldlayout->addWidget(emptymacroscopicfieldLabel,1,0,Qt::AlignLeft);
  studynearfieldlayout->addWidget(macroscopicfieldLabel,1,1,Qt::AlignLeft);
  studynearfieldlayout->addWidget(macroscopicfield,1,2,Qt::AlignLeft);
  emptyrangeofstudyLabel = new QLabel(" ");
  rangeofstudyLabel = new QLabel("Range of study");
  rangeofstudy = new QComboBox();
  rangeofstudy->addItems(options->rangeofstudyList);
  rangeofstudy->setCurrentIndex(options->getNproche());
  studynearfieldlayout->addWidget(emptyrangeofstudyLabel,2,0,Qt::AlignLeft);
  studynearfieldlayout->addWidget(rangeofstudyLabel,2,1,Qt::AlignLeft);
  studynearfieldlayout->addWidget(rangeofstudy,2,2,Qt::AlignLeft);

  nxm = new QLineEdit(QString::number(options->getNxm()));
  nxm->setFixedWidth(40);
  nym = new QLineEdit(QString::number(options->getNym()));
  nym->setFixedWidth(40);
  nzm = new QLineEdit(QString::number(options->getNzm()));
  nzm->setFixedWidth(40);

  QBoxLayout   *beamlayout = new QBoxLayout(QBoxLayout::LeftToRight);
  beamlayout->addWidget(beam);
  beamlayout->addWidget(beamButton);

  QBoxLayout   *objectlayout = new QBoxLayout(QBoxLayout::LeftToRight);
  objectlayout->addWidget(object);
  objectlayout->addWidget(objectButton);

  QBoxLayout   *epsilonlayout = new QBoxLayout(QBoxLayout::LeftToRight);
  epsilonlayout->addWidget(anisotropy);
  epsilonlayout->addWidget(epsilonButton);

  QFrame *hsep000 = new QFrame(this);
  hsep000->setFrameShape(QFrame::HLine);
  hsep000->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep000);
  layout->addRow(nrigLabel->text(),nrig);
  layout->addRow(nreadLabel->text(),nread);
  layout->addRow(startnreadlayout);
  layout->addRow(nmatlabLabel->text(),nmatlab);
  layout->addRow(calculationlayout);
  QFrame *hsep0 = new QFrame(this);
  QFrame *hsep00 = new QFrame(this);
  hsep0->setFrameShape(QFrame::HLine);
  hsep0->setFrameShadow(QFrame::Sunken);
  hsep00->setFrameShape(QFrame::HLine);
  hsep00->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep0);
  layout->addRow("Illumination properties",new QLabel(""));
  layout->addRow(hsep00);
  layout->addRow(wavelengthLabel,wavelength);
  layout->addRow(P0Label,P0);
  layout->addRow(W0Label,W0);
  layout->addRow(beamLabel,beamlayout);
  layout->addRow(wavemultinumberLabel,wavemultinumber);
  QFrame *hsep01 = new QFrame(this);
  QFrame *hsep001 = new QFrame(this);
  hsep01->setFrameShape(QFrame::HLine);
  hsep01->setFrameShadow(QFrame::Sunken);
  hsep001->setFrameShape(QFrame::HLine);
  hsep001->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep01);
  layout->addRow("Object properties",new QLabel(""));
  layout->addRow(hsep001);
  layout->addRow(objectLabel,objectlayout);
  layout->addRow(objectnumberLabel,objectnumber);
  layout->addRow(anisotropyLabel,epsilonlayout);
  layout->addRow(discretizationLabel,discretization);
  QFrame *hsep02 = new QFrame(this);
  QFrame *hsep002 = new QFrame(this);
  hsep02->setFrameShape(QFrame::HLine);
  hsep02->setFrameShadow(QFrame::Sunken);
  hsep002->setFrameShape(QFrame::HLine);
  hsep002->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep02);
  layout->addRow("Study",new QLabel(""));
  layout->addRow(hsep002);
  dipolepsilon = new QCheckBox(this);
  connect(dipolepsilon, SIGNAL(stateChanged(int)),this,
	SLOT(dipolepsilonCheckBoxStateChanged(int)));
  layout->addRow("Only dipoles with epsilon",dipolepsilon);
  farfield = new QCheckBox(this);
  connect(farfield , SIGNAL(stateChanged(int)),this,
	SLOT(farfieldCheckBoxStateChanged(int)));
  layout->addRow("Far field",farfield);
  layout->addRow(studyfarfieldlayout);
  force = new QCheckBox(this);
  connect(force , SIGNAL(stateChanged(int)),this,
	SLOT(forceCheckBoxStateChanged(int)));
  layout->addRow("Force",force);
  layout->addRow(studyforcelayout);
  nearfield = new QCheckBox(this);
  connect(nearfield , SIGNAL(stateChanged(int)),this,
	SLOT(nearfieldCheckBoxStateChanged(int)));
  layout->addRow("Near field",nearfield);
  layout->addRow(studynearfieldlayout);
  QFrame *hsep1 = new QFrame(this);
  QFrame *hsep2 = new QFrame(this);
  hsep1->setFrameShape(QFrame::HLine);
  hsep1->setFrameShadow(QFrame::Sunken);
  hsep2->setFrameShape(QFrame::HLine);
  hsep2->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep1);
  layout->addRow("Numerical parameters",new QLabel(""));
  layout->addRow(hsep2);
  layout->addRow(toleranceLabel,tolerance);
  layout->addRow(methodeitLabel,methodeit);
  layout->addRow("Nxm:",nxm);
  layout->addRow("Nym:",nym);
  layout->addRow("Nzm:",nzm);
  layout->addRow(polarizabilityLabel,polarizability);
  layout->addRow(quadLabel,quad);
  layout->addRow(nfft2dLabel,nfft2d);
  this->setLayout(layout);
  update();
  beamconfigdlg = NULL;
  objectconfigdlg = NULL;
  epsilonconfigdlg = NULL;

  // init now object combo index (triggers the handle)
  QLOG_INFO ( ) << "Number of objects : " << options->getObjectNumber();
  objectnumber->setValue(options->getObjectNumber());
  object->setCurrentIndex(object->findText(options->getObject()));
  connect(objectnumber , SIGNAL(valueChanged(int)),this,
	SLOT(handleObjectNumberSelectionChanged(int)));

  // init now wavemulti combo index (triggers the handle)
  QLOG_INFO ( ) << "Number of plane waves : " << options->getWaveMultiNumber();
  wavemultinumber->setValue(options->getWaveMultiNumber());
  beam->setCurrentIndex(beam->findText(options->getBeam()));
  connect(wavemultinumber , SIGNAL(valueChanged(int)),this,
	SLOT(handleWaveMultiNumberSelectionChanged(int)));
  
}
OptionsWidget::~OptionsWidget()
{
  QLOG_DEBUG ( ) << "Deleting OptionsWidget"; 
}
void 
OptionsWidget::execute() {
  RunWidget *runwidget = NULL;
  OptionsWidget *optionswidget = NULL;
  if (mainwindow) {
    QTabWidget *tabwidget = (QTabWidget*) mainwindow->findChild<QTabWidget*>("TabWidget");
    QMainWindow *currentWindow = (QMainWindow*)tabwidget->currentWidget();
    runwidget = currentWindow->findChild<RunWidget *>("Run");
    optionswidget = currentWindow->findChild<OptionsWidget *>("Options");
    optionswidget->updateOptions();
    runwidget->execute();
   }
}
void
OptionsWidget::saveName() {
  QLOG_INFO() << "PRESENT SAVE NAME " << name;
  this->updateOptions();
  if (name == "new") {
    nameLineEdit->clear();
    descriptionLineEdit->clear();
    cfgWindow->exec();
    QLOG_INFO() << "NEW SAVE NAME " << name;
  }
  else {
   nameLineEdit->setText(name);
   this->finish();
  }
}
void 
OptionsWidget::finish() {
 name = nameLineEdit->text();
 QLOG_INFO() << "OptionsWidget::finish() " << name;
 if ( name == "" ) {
   QMessageBox::warning(this, "Warning:", "Enter a valid name");
   return;
 }
 options->setName(name);
 options->setDescription(descriptionLineEdit->text());
 description = options->getDescription();
 options->saveDb(name,description);
 cfgWindow->hide();
 if (mainwindow) {
   QTabWidget *tabwidget = (QTabWidget*) mainwindow->findChild<QTabWidget*>("TabWidget");
   QMainWindow *currentWindow = (QMainWindow*)tabwidget->currentWidget();
   QDockWidget *dockWidget = currentWindow->findChild<QDockWidget*>("DockWidget");
   RunWidget *runwidget = currentWindow->findChild<RunWidget *>("Run");
   QDockWidget *dockWidgetCentral = runwidget->findChild<QDockWidget*>("DockWidgetCentral");
   QDockWidget *dockWidgetOutput = runwidget->findChild<QDockWidget*>("DockWidgetOutput");
   QDockWidget *dockWidgetPlot = runwidget->findChild<QDockWidget*>("DockWidgetPlot");
   dockWidget->setWindowTitle(name);
   dockWidgetCentral->setWindowTitle(name);
   dockWidgetOutput->setWindowTitle(name);
   dockWidgetPlot->setWindowTitle(name);
   tabwidget->setTabText(tabwidget->currentIndex(),name);
 }
  emit updateOptionsWindow();
}
void 
OptionsWidget::cancel() {
 cfgWindow->hide();
}
void 
OptionsWidget::handleObjectSelectionChanged(int index){

 QLOG_INFO() << "OptionsWidget::handleObjectSelectionChanged";
 if ( object->currentText() != "nspheres" && 
      object->currentText() != "concentricsphere")  {
   this->setObjectNumber(1);
   objectnumber->setEnabled( false );
 }
 else {
   objectnumber->setValue(options->getObjectNumber());
   objectnumber->setEnabled( true );
 }
}
void 
OptionsWidget::handleObjectNumberSelectionChanged(int index) {
//if ( object->currentText() == "nspheres" ||
//     object->currentText() == "concentricsphere")
   options->setObjectNumber(this->getObjectNumber());
}
void 
OptionsWidget::handleBeamSelectionChanged(int index){

 QLOG_INFO() << "OptionsWidget::handleBeamSelectionChanged";
 if ( beam->currentText() != "wavelinearmulti" )  {
   this->setWaveMultiNumber(1);
   wavemultinumber->setEnabled( false );
 }
 else {
   wavemultinumber->setValue(options->getWaveMultiNumber());
   wavemultinumber->setEnabled( true );
 }
}
void 
OptionsWidget::handleWaveMultiNumberSelectionChanged(int index) {
//if ( object->currentText() == "nspheres" ||
//     object->currentText() == "concentricsphere")
   options->setWaveMultiNumber(this->getWaveMultiNumber());
}
void 
OptionsWidget::nreadCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   filerereadLabel->show();
   filereread->show();
  }
  else if (state == Qt::Unchecked) {
   filerereadLabel->hide();
   filereread->hide(); 
  }
}
void
OptionsWidget::nmatlabCheckBoxStateChanged(int state) {
}
void 
OptionsWidget::dipolepsilonCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   farfield->setChecked(false);
   nearfield->setChecked(false);
   force->setChecked(false);
  }
}
void 
OptionsWidget::farfieldCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   emptycrosssectionLabel->show();
   crosssectionLabel->show();
   crosssection->show();
   emptycrosssectionpoyntingLabel->show();
   crosssectionpoyntingLabel->show();
   crosssectionpoynting->show();
   emptynenergieLabel->show();
   nenergieLabel->show();
   nenergie->show();
   emptymicroscopyLabel->show();  
   microscopyLabel->show();
   microscopy->show();
   dipolepsilon->setChecked(false);
  }
  else if (state == Qt::Unchecked) {
   emptycrosssectionLabel->hide();
   crosssection->setChecked(false);
   crosssectionLabel->hide();
   crosssection->hide();  
   emptycrosssectionpoyntingLabel->hide();
   crosssectionpoynting->setChecked(false);
   crosssectionpoyntingLabel->hide();
   crosssectionpoynting->hide();
   crosssectionpoyntingCheckBoxStateChanged(Qt::Unchecked);
   emptynenergieLabel->hide(); 
   nenergie->setChecked(false);
   nenergieLabel->hide();
   nenergie->hide();
   nenergieCheckBoxStateChanged(Qt::Unchecked);
   emptymicroscopyLabel->hide();  
   microscopy->setChecked(false);
   microscopyLabel->hide();
   microscopy->hide();
   microscopyCheckBoxStateChanged(Qt::Unchecked);
  }
}
void 
OptionsWidget::crosssectionCheckBoxStateChanged(int state) {
}
void 
OptionsWidget::nenergieCheckBoxStateChanged(int state) {
}
void 
OptionsWidget::crosssectionpoyntingCheckBoxStateChanged(int state) {
QLOG_INFO() << "OptionsWidget::crosssectionpoyntingCheckBoxStateChanged" << state;
if (state == Qt::Checked) {
   quickdiffractLabel->show();
   quickdiffract->show();
   nthetaLabel->show();
   ntheta->show();
   nphiLabel->show();
   nphi->show();
  }
  else if (state == Qt::Unchecked) {
   quickdiffractLabel->hide();
   quickdiffract->setChecked(false);
   quickdiffract->hide();
   nthetaLabel->hide();
   ntheta->hide();
   nphiLabel->hide();
   nphi->hide();
  }
}
void 
OptionsWidget::microscopyCheckBoxStateChanged(int state) {
if (state == Qt::Checked) {
   microscopyFFTLabel->show();
   microscopyFFT->show();
   naLabel->show();
   na->show();
   grossLabel->show();
   gross->show();
  }
  else if (state == Qt::Unchecked) {
   microscopyFFTLabel->hide();
   microscopyFFT->setChecked(false);
   microscopyFFT->hide();
   naLabel->hide();
   na->hide();
   grossLabel->hide();
   gross->hide();
  }
}
void 
OptionsWidget::forceCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   emptyopticalforceLabel->show();
   opticalforceLabel->show();
   opticalforce->show();
   emptyopticalforcedensityLabel->show();
   opticalforcedensityLabel->show();
   opticalforcedensity->show();
   emptyopticaltorqueLabel->show();
   opticaltorqueLabel->show();
   opticaltorque->show();
   emptyopticaltorquedensityLabel->hide();
   opticaltorquedensityLabel->show();
   opticaltorquedensity->show();
   dipolepsilon->setChecked(false);
  }
  else if (state == Qt::Unchecked) {
   emptyopticalforceLabel->hide();
   opticalforce->setChecked(false);
   opticalforceLabel->hide();
   opticalforce->hide();
   emptyopticalforcedensityLabel->hide();
   opticalforcedensity->setChecked(false);
   opticalforcedensityLabel->hide();
   opticalforcedensity->hide();
   emptyopticaltorqueLabel->hide();
   opticaltorque->setChecked(false);
   opticaltorqueLabel->hide();
   opticaltorque->hide();
   emptyopticaltorquedensityLabel->hide();
   opticaltorquedensity->setChecked(false);
   opticaltorquedensityLabel->hide();
   opticaltorquedensity->hide();
  }
}
void 
OptionsWidget::opticalforceCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked) {
    opticalforcedensity->setChecked(false);
    opticaltorque->setChecked(false);
    if ( opticaltorque->isChecked() == false &&
        opticaltorquedensity->isChecked() == false )
     force->setChecked(false);
   }
}
void 
OptionsWidget::opticalforcedensityCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked) {
    if (opticalforce->isChecked() == false &&
        opticaltorque->isChecked() == false &&
        opticaltorquedensity->isChecked() == false )
     force->setChecked(false);
   }
   else if (state == Qt::Checked)
     opticalforce->setChecked(true);
   
}
void 
OptionsWidget::opticaltorqueCheckBoxStateChanged(int state) {
    if (state == Qt::Unchecked) {
      opticaltorquedensity->setChecked(false);
      if ( opticalforcedensity->isChecked() == false &&
        opticalforce->isChecked() == false )
     force->setChecked(false);
   }
   else if (state == Qt::Checked)
     opticalforce->setChecked(true);
}
void 
OptionsWidget::opticaltorquedensityCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked) {
    if (opticalforcedensity->isChecked() == false &&
        opticaltorque->isChecked() == false &&
        opticalforce->isChecked() == false )
      force->setChecked(false);
   }
   else if (state == Qt::Checked)
     opticaltorque->setChecked(true);
}
void 
OptionsWidget::nearfieldCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   emptyrangeofstudyLabel->show();
   rangeofstudyLabel->show();
   rangeofstudy->show();
   emptylocalfieldLabel->show();
   localfieldLabel->show();
   localfield->show();
   emptymacroscopicfieldLabel->show();
   macroscopicfieldLabel->show();
   macroscopicfield->show();
   dipolepsilon->setChecked(false);
   nxm->setText(QString::number(this->getDiscretization()));
   nym->setText(QString::number(this->getDiscretization()));
   nzm->setText(QString::number(this->getDiscretization()));
   nxm->show();
   nym->show();
   nzm->show();
  }
  else if (state == Qt::Unchecked) {
   emptylocalfieldLabel->hide();
   localfield->setChecked(false);
   localfield->hide();  
   localfieldLabel->hide();
   emptymacroscopicfieldLabel->hide();
   macroscopicfield->setChecked(false);
   macroscopicfield->hide();
   macroscopicfieldLabel->hide();
   emptyrangeofstudyLabel->hide();
   rangeofstudy->hide();  
   rangeofstudyLabel->hide();
   nxm->setText(QString::number(this->getDiscretization()));
   nym->setText(QString::number(this->getDiscretization()));
   nzm->setText(QString::number(this->getDiscretization()));
   nxm->hide();
   nym->hide();
   nzm->hide();
  }
}
void 
OptionsWidget::localfieldCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked)
    if (macroscopicfield->isChecked() == false )
     nearfield->setChecked(false);
}
void 
OptionsWidget::macroscopicfieldCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked)
    if (localfield->isChecked() == false )
     nearfield->setChecked(false);
}
QString 
OptionsWidget::getFilereread(){
  return filereread->text();
}
double 
OptionsWidget::getWavelength(){
  return wavelength->text().toDouble();
}
double 
OptionsWidget::getP0(){
  return P0->text().toDouble();
}
double 
OptionsWidget::getW0(){
  return W0->text().toDouble();
}
QString 
OptionsWidget::getBeam(){
  return beam->currentText();
}
QString 
OptionsWidget::getObject(){
  return object->currentText();
}
int 
OptionsWidget::getObjectNumber(){
  return objectnumber->value();
}
int 
OptionsWidget::getWaveMultiNumber(){
  return wavemultinumber->value();
}
QString 
OptionsWidget::getAnisotropy(){
  return anisotropy->currentText();
}

int 
OptionsWidget::getDiscretization(){
  return discretization->text().toInt();
}
double 
OptionsWidget::getTolerance(){
  return tolerance->text().toDouble();
}
QString 
OptionsWidget::getMethodeit(){
  return methodeit->currentText();
}
QString 
OptionsWidget::getPolarizability(){
  return polarizability->currentText();
}
int 
OptionsWidget::getQuad(){
  return quad->currentText().toInt();
}
int 
OptionsWidget::getnfft2d(){
  return nfft2d->currentText().toInt();
}
void 
OptionsWidget::setWavelength(double _wavelength){
  wavelength->setText(QString::number(_wavelength));
}
void 
OptionsWidget::setP0(double _P0){
  P0->setText(QString::number(_P0));
}
void 
OptionsWidget::setW0(double _W0){
  W0->setText(QString::number(_W0));
}
void 
OptionsWidget::setBeam(QString _beam){
  beam->setCurrentIndex(beam->findText(_beam));
}
void 
OptionsWidget::setObject(QString _object){
  object->setCurrentIndex(object->findText(_object));
}
void 
OptionsWidget::setObjectNumber(int _objectnumber){
  objectnumber->setValue(_objectnumber);
}
void 
OptionsWidget::setWaveMultiNumber(int _wavemultinumber){
  wavemultinumber->setValue(_wavemultinumber);
}
void 
OptionsWidget::setAnisotropy(QString _anisotropy){
anisotropy->setCurrentIndex(anisotropy->findText(_anisotropy));
}
void 
OptionsWidget::setDiscretization(int _discretization){
discretization->setText(QString::number(_discretization));
}
void 
OptionsWidget::setTolerance(double _tolerance){
tolerance->setText(QString::number(_tolerance));
}
void 
OptionsWidget::setMethodeit(QString _methodeit){
methodeit->setCurrentIndex(methodeit->findText(_methodeit));
}
void 
OptionsWidget::setPolarizability(QString _polarizability){
polarizability->setCurrentIndex(polarizability->findText(_polarizability));
}
void 
OptionsWidget::setQuad(int _quad){
quad->setCurrentIndex(quad->findText(QString::number(_quad)));
}
void 
OptionsWidget::setnfft2d(int _nfft2d){
nfft2d->setCurrentIndex(nfft2d->findText(QString::number(_nfft2d)));
}
void OptionsWidget::configureBeam() {
   this->updateOptions();
   if (beamconfigdlg) delete beamconfigdlg;
   beamconfigdlg = new BeamConfigDialog(this,Qt::Widget,options);
   beamconfigdlg->exec();
}
void OptionsWidget::configureObject(){
   this->updateOptions();
   if (objectconfigdlg) delete objectconfigdlg;
   objectconfigdlg = new ObjectConfigDialog(this,Qt::Widget,options);
   objectconfigdlg->exec();
}
void OptionsWidget::configureEpsilon() {
   this->updateOptions();
   if (epsilonconfigdlg) delete epsilonconfigdlg;
   epsilonconfigdlg = new EpsilonConfigDialog(this,Qt::Widget,options);
   epsilonconfigdlg->exec();
}
void 
OptionsWidget::updateFarfield() {
  crosssection->setChecked(options->getCrosssection());
  crosssectionpoynting->setChecked(options->getCrosssectionpoynting());
  quickdiffract->setChecked(options->getQuickdiffract());
  nenergie->setChecked(options->getNenergie());
  microscopy->setChecked(options->getMicroscopy());
  microscopyFFT->setChecked(options->getMicroscopyFFT());
  ntheta->setText(QString::number(options->getNtheta()));
  nphi->setText(QString::number(options->getNphi()));
  na->setText(QString::number(options->getNA()));
  gross->setText(QString::number(options->getGross()));
  if (crosssection->isChecked() || crosssectionpoynting->isChecked() ||
      microscopy->isChecked() || nenergie->isChecked())
   farfield->setChecked(true);
  farfieldCheckBoxStateChanged(farfield->isChecked());
}
void 
OptionsWidget::updateNearfield() {
  localfield->setChecked(options->getLocalfield());
  macroscopicfield->setChecked(options->getMacroscopicfield());
  rangeofstudy->setCurrentIndex(options->getNproche());
  if (localfield->isChecked() || macroscopicfield->isChecked())
   nearfield->setChecked(true);
  nearfieldCheckBoxStateChanged(nearfield->isChecked()); 
}
void 
OptionsWidget::updateForce() {
  opticalforce->setChecked(options->getOpticalforce());
  opticalforcedensity->setChecked(options->getOpticalforcedensity());
  opticaltorque->setChecked(options->getOpticaltorque());
  opticaltorquedensity->setChecked(options->getOpticaltorquedensity());
  if (opticalforce->isChecked() || opticalforcedensity->isChecked() ||
      opticaltorque->isChecked() || opticaltorquedensity->isChecked())
   force->setChecked(true);
  forceCheckBoxStateChanged(force->isChecked());
}
void 
OptionsWidget::update() {
  filereread->setText(options->getFilereread());
  this->setWavelength(options->getWavelength());
  this->setP0(options->getP0());
  this->setW0(options->getW0());
  this->setBeam(options->getBeam());
  this->setObject(options->getObject());
  this->setObjectNumber(options->getObjectNumber());
  this->setWaveMultiNumber(options->getWaveMultiNumber());
  this->setAnisotropy(options->getAnisotropy());
  this->setDiscretization(options->getDiscretization());
  this->setTolerance(options->getTolerance());
  this->setMethodeit(options->getMethodeit());
  this->setPolarizability(options->getPolarizability());
  this->setQuad(options->getQuad());
  this->setnfft2d(options->getnfft2d());
  QLOG_DEBUG() << " NREAD " << options->getNread();
  nread->setChecked(options->getNread());
  nreadCheckBoxStateChanged(nread->isChecked());
  nmatlab->setChecked(options->getNmatlab());
  nmatlabCheckBoxStateChanged(nmatlab->isChecked());
  dipolepsilon->setChecked(options->getDipolepsilon());
  updateFarfield();
  updateForce();
  updateNearfield();
  nxm->setText(QString::number(options->getNxm()));
  nym->setText(QString::number(options->getNym()));
  nzm->setText(QString::number(options->getNzm()));
}
void 
OptionsWidget::updateOptions() {
  options->setNread(nread->isChecked());
  options->setNmatlab(nmatlab->isChecked());
  options->setWavelength(this->getWavelength());
  options->setP0(this->getP0());
  options->setW0(this->getW0());
  options->setBeam(this->getBeam());
  options->setObject(this->getObject());
  options->setObjectNumber(this->getObjectNumber());
  options->setWaveMultiNumber(this->getWaveMultiNumber());
  options->setAnisotropy(this->getAnisotropy());
  options->setDiscretization(this->getDiscretization());
  options->setTolerance(this->getTolerance());
  options->setMethodeit(this->getMethodeit());
  options->setPolarizability(this->getPolarizability());
  options->setQuad(this->getQuad());
  options->setDipolepsilon(dipolepsilon->isChecked());
  options->setFarfield(farfield->isChecked());
  options->setNearfield(nearfield->isChecked());
  options->setForce(force->isChecked());
  options->setLocalfield(localfield->isChecked());
  options->setMacroscopicfield(macroscopicfield->isChecked());
  options->setCrosssection(crosssection->isChecked());
  options->setCrosssectionpoynting(crosssectionpoynting->isChecked());
  options->setQuickdiffract(quickdiffract->isChecked());
  options->setNrig(nrig->currentIndex());
  options->setNenergie(nenergie->isChecked());
  options->setMicroscopy(microscopy->isChecked());
  options->setMicroscopyFFT(microscopyFFT->isChecked());
  options->setNenergie(nenergie->isChecked());
  options->setOpticalforce(opticalforce->isChecked());
  options->setOpticalforcedensity(opticalforcedensity->isChecked());
  options->setOpticaltorque(opticaltorque->isChecked());
  options->setOpticaltorquedensity(opticaltorquedensity->isChecked());
  options->setNproche(rangeofstudy->currentIndex());
  QLOG_DEBUG() << "NPROCHE " << options->getNproche();

  if ( options->getNread() == true )
   options->setFilereread(this->getFilereread());
  if ( options->getNearfield() == true ) {
    options->setNxm(nxm->text().toInt());
    options->setNym(nym->text().toInt());
    options->setNzm(nzm->text().toInt());
  }
  else {
    QLOG_INFO() << " UPDATE NXM = " << this->getDiscretization();
    options->setNxm(this->getDiscretization());
    options->setNym(this->getDiscretization());
    options->setNzm(this->getDiscretization());
  }
  options->setNtheta(ntheta->text().toInt());
  options->setNphi(nphi->text().toInt());
  options->setNA(na->text().toDouble());
  options->setGross(gross->text().toDouble());
  options->setnfft2d(this->getnfft2d());  
  if ( options->getLocalfield() == false && options->getMacroscopicfield() == false)
    options->setNearfield(false);
  if ( options->getOpticalforce() == false && options->getOpticalforcedensity() == false &&
        options->getOpticaltorque() == false && options->getOpticaltorquedensity() == false )
    options->setForce(false);
  if ( options->getNenergie() == false && 
       options->getMicroscopy() == false && 
       options->getCrosssection() == false &&
       options->getCrosssectionpoynting() == false)
    options->setFarfield(false);
  QLOG_DEBUG() << "OptionsWidget::updateOptions> " << QString::number(options->getNA());
  QLOG_DEBUG() << "OptionsWidget::updateOptions> " << QString::number(options->getGross());
  QLOG_DEBUG() << "OptionsWidget::updateOptions> " << QString::number(options->getNrig());
  QLOG_DEBUG() << "OptionsWidget::updateOptions> crosssection " << crosssection->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> crosssectionpoynting " << crosssectionpoynting->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> quickdiffract " << quickdiffract->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> localfield " << localfield->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> macroscopicfield " << macroscopicfield->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticalforce " << opticalforce->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticalforcedensity " << opticalforcedensity->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticaltorque " << opticaltorque->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticaltorquedensity " << opticaltorquedensity->isChecked();
  QLOG_INFO() << "OptionsWidget::updateOptions> nenergie " << nenergie->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> microscopy " << microscopy->isChecked();
}

