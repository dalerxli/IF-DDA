#include "cdmOptionsWindow.h"
#include "cdmMain.h"

OptionsWindow::OptionsWindow( QMainWindow *_parent, Options *_options, Qt::WindowFlags fl)
  : QMainWindow( _parent, fl )
{
  parent = _parent;
  optionswidget = NULL;
  optionstable = NULL;
  options = _options;
  
  QWidget *centralWidget = new QWidget(this);
  this->setCentralWidget(centralWidget);
  vboxlayout = new QVBoxLayout();
  centralWidget->setLayout(vboxlayout);
  this->setWindowTitle("Configurations saved");
  this->setMinimumHeight(300);
  this->setMinimumWidth(400);
  
  optionstitle = new QLabel("Select a configuration:");
  vboxlayout->addWidget(optionstitle);
  optionsview = new QTableView();
  QPalette palette = optionsview->palette();
  palette.setBrush(QPalette::Base, Qt::transparent);
  optionsview->setPalette(palette);
  optionsview->setItemDelegate(new ColorDelegate());
  optionsview->setEditTriggers(QAbstractItemView::AllEditTriggers);
  //optionsview->setSelectionBehavior(QAbstractItemView::SelectRows);
  vboxlayout->addWidget(optionsview);

  saveButton = new QPushButton("Save",this);
  saveButton->setToolTip("Enter a new name or select a line to save the current calculation");
  saveButton->setFixedHeight(30);
  saveButton->setFixedWidth(150);
  QObject::connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  vboxlayout->addWidget(saveButton);
  removeButton = new QPushButton("Delete",this);
  removeButton->setToolTip("Select a line to delete");
  removeButton->setFixedHeight(30);
  removeButton->setFixedWidth(150);
  QObject::connect(removeButton, SIGNAL(clicked()), this, SLOT(remove()));
  vboxlayout->addWidget(removeButton);
  loadButton = new QPushButton("Load",this);
  loadButton->setToolTip("Select a line to load");
  loadButton->setFixedHeight(30);
  loadButton->setFixedWidth(150);
  QObject::connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
  vboxlayout->addWidget(loadButton);
  exportButton = new QPushButton("Export",this);
  exportButton->setToolTip("Select a line to export");
  exportButton->setFixedHeight(30);
  exportButton->setFixedWidth(150);
  QObject::connect(exportButton, SIGNAL(clicked()), this, SLOT(tofile()));
  vboxlayout->addWidget(exportButton);
  initOptionsTbl();
} 

OptionsWindow::~OptionsWindow()
{
  QLOG_DEBUG ( ) << "Deleting OptionsWindow";
  if (optionstable){ delete optionstable; optionstable = NULL;}
}

void OptionsWindow::closeEvent(QCloseEvent* event)
{
  event->accept();  
  QLOG_DEBUG ( ) << "Closing OptionsWindow";
  this->hide();
}
void OptionsWindow::initOptionsTbl() {
  QTabWidget *tabwidget = NULL;
  QMainWindow *currentWindow = NULL;
  tabwidget = (QTabWidget*) parent->findChild<QTabWidget*>("TabWidget");
  if (tabwidget) {
   currentWindow = (QMainWindow*)tabwidget->currentWidget();
   optionswidget = currentWindow->findChild<OptionsWidget *>("Options");
  }
  /*if (optionswidget)
    options = optionswidget->options;
  else {
   if (options) {delete options; options = NULL;}
   options = new Options();
   options->initDb();
  } */
  if (optionstable) {delete optionstable; optionstable = NULL;}
  dbpath = "options.db3";
  optionstable = new ManifestModel(this,QSqlDatabase::database(dbpath));
  optionstable->setTable("options_tbl");
  optionstable->setEditStrategy(QSqlTableModel::OnManualSubmit);
  optionsview->setModel(optionstable);
  optionstable->setFilter("name not like 'new'");
  optionstable->select();
  optionstable->insertRow(optionstable->rowCount());
  optionsview->resizeColumnsToContents();
  optionsview->resizeRowsToContents();
  optionsview->setCurrentIndex (optionsview->model()->index(optionstable->rowCount() - 1, 
		    0, QModelIndex()));
}
void 
OptionsWindow::save(){
  if (optionsview->selectionModel()->hasSelection() == false) {
     this->showWarning("Select a valid configuration to save");
     return;
  }
  QModelIndex selected = optionsview->selectionModel()->currentIndex();
  QSqlRecord record = optionstable->record(selected.row());
  QString name = record.value(0).toString();
  QString description = record.value(1).toString();
  if (name =="") {
   this->showWarning("Select a valid configuration to save");
   return;
  }
  QTabWidget *tabwidget = (QTabWidget*) parent->findChild<QTabWidget*>("TabWidget");
  QMainWindow *currentWindow = (QMainWindow*)tabwidget->currentWidget();
  optionswidget = currentWindow->findChild<OptionsWidget *>("Options");
  if (optionswidget)
    optionswidget->updateOptions();
  else
    options->loadDb(name);
  options->saveDb(name,description);
  initOptionsTbl();
}
void 
OptionsWindow::load(){
  if (optionsview->selectionModel()->hasSelection() == false) {
    this->showWarning("Select a valid configuration to load");
  return;
  }
  QModelIndexList selectedlist = optionsview->selectionModel()->selectedIndexes();
  for (int i = selectedlist.size() - 1 ; i >=0 ; i--) {
    QModelIndex selected = selectedlist.at(i);
    QSqlRecord record = optionstable->record(selected.row());
    QString name = record.value(0).toString();
    options->loadDb(name);
    CdmMain *mainwindow = (CdmMain*) parent; 
    mainwindow->openWidget(options);
 }
 this->hide();
}
void 
OptionsWindow::remove(){
  if (optionsview->selectionModel()->hasSelection() == false) {
    this->showWarning("Select a valid configuration to delete");
    return;
  }
  QModelIndexList selectedlist = optionsview->selectionModel()->selectedIndexes();
  for (int i = selectedlist.size() - 1 ; i >=0 ; i--) {
    QModelIndex selected = selectedlist.at(i);
    QSqlRecord record = optionstable->record(selected.row());
    QString name = record.value(0).toString();
    QString description = record.value(1).toString();
    if (name =="") {
      this->showWarning("Select a valid configuration to delete");
      return;
    }
    options->removeDb(name);
  }
  initOptionsTbl();
}
void 
OptionsWindow::updatewindow(){
  initOptionsTbl();
}
void 
OptionsWindow::tofile(){
 if (optionsview->selectionModel()->hasSelection() == false) {
    this->showWarning("Select a valid configuration to export");
    return;
  }
  QModelIndex selected = optionsview->selectionModel()->currentIndex();
  QSqlRecord record = optionstable->record(selected.row());
  QString name = record.value(0).toString();
  if (name =="") {
   this->showWarning("Select a valid configuration to export");
   return;
  }
  options->loadDb(name);
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                     name + ".opt", tr("Options (*.opt)"));
  QFile optfile(fileName);
  if (!optfile.open(QIODevice::WriteOnly | QIODevice::Text))
     return;
  QTextStream opt(&optfile);
  // Fill ASCII options here
  opt << "calculation options [0=rigorous, 1=renormalized Born, 2=Born, 3=Born order, 4=renormalized Rytov, 5=Rytov, 6=Beam propagation, 7=renormalized Beam propagation]:" << options->getNrig() << endl;
  opt << "illumination properties options:" << endl;
  opt << " wavelength:" << options->getWavelength() << endl;
  opt << " P0:" << options->getP0() << endl;
  opt << " W0:" << options->getW0() << endl;
  opt << " beam:" << options->getBeam() << endl;
  if (options->getBeam() == "pwavecircular") {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization L (-1) R (1):" << options->getPolarizationRL() << endl;
  }
  else if (options->getBeam() == "pwavelinear") {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization TM (1) TE (0):" << options->getPolarizationTM() << endl;
  }
  else if (options->getBeam() == "gwavecircular" || options->getBeam() == "gparawavecircular" ) {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization L (-1) R (1):" << options->getPolarizationRL() << endl;
    opt << "  gaussian X:" << options->getXgaus() << endl;
    opt << "  gaussian Y:" << options->getYgaus() << endl;
    opt << "  gaussian Z:" << options->getZgaus() << endl;
  }
  else if (options->getBeam() == "gwavelinear" || options->getBeam() == "gparawavelinear" ) {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization TM (1) TE (0):" << options->getPolarizationTM() << endl;
    opt << "  gaussian X:" << options->getXgaus() << endl;
    opt << "  gaussian Y:" << options->getYgaus() << endl;
    opt << "  gaussian Z:" << options->getZgaus() << endl;
  }
  
  opt << "object properties options:" << endl;
  for (int i = 0 ; i < options->getObjectNumber(); i++) {
    opt << " object" << i << ":" << options->getObject() << endl;
    if (options->getObject() == "sphere" || options->getObject() == "nspheres") {
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
    }
    else if (options->getObject() == "inhomosphere") {
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  coherence length:" << options->getSpherecoherencelength() << endl;
      opt << "  standard deviation:" << options->getSpherestandardev() << endl;
    }
     else if (options->getObject() == "randomsphere") {
      opt << " cube side X:" << options->getCubesidex() << endl;
      opt << " cube side Y:" << options->getCubesidey() << endl;
      opt << " cube side Z:" << options->getCubesidez() << endl;
      opt << " position X:" << options->getPositionx().at(i) << endl;
      opt << " position Y:" << options->getPositiony().at(i) << endl;
      opt << " position Z:" << options->getPositionz().at(i) << endl;
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  density:" << options->getDensity() << endl;
     }
    else if (options->getObject() == "concentricsphere") {
      opt << " radius:" << options->getSphereradius().at(i) << endl;
      if (i == 0) {
        opt << " position X:" << options->getPositionx().at(i) << endl;
        opt << " position Y:" << options->getPositiony().at(i) << endl;
        opt << " position Z:" << options->getPositionz().at(i) << endl;
      }
    }
    else if (options->getObject() == "cube") {
      opt << " cube side:" << options->getCubeside() << endl;
      opt << " position X:" << options->getPositionx().at(i) << endl;
      opt << " position Y:" << options->getPositiony().at(i) << endl;
      opt << " position Z:" << options->getPositionz().at(i) << endl;
    }
    else if (options->getObject() == "inhomocuboid") {
      opt << " cube side X:" << options->getCubesidex() << endl;
      opt << " cube side Y:" << options->getCubesidey() << endl;
      opt << " cube side Z:" << options->getCubesidez() << endl;
      opt << " position X:" << options->getPositionx().at(i) << endl;
      opt << " position Y:" << options->getPositiony().at(i) << endl;
      opt << " position Z:" << options->getPositionz().at(i) << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  coherence length:" << options->getSpherecoherencelength() << endl;
      opt << "  standard deviation:" << options->getSpherestandardev() << endl;
    }
      else if (options->getObject() == "cuboid") {
      opt << " cube side X:" << options->getCubesidex() << endl;
      opt << " cube side Y:" << options->getCubesidey() << endl;
      opt << " cube side Z:" << options->getCubesidez() << endl;
      opt << " position X:" << options->getPositionx().at(i) << endl;
      opt << " position Y:" << options->getPositiony().at(i) << endl;
      opt << " position Z:" << options->getPositionz().at(i) << endl;
      opt << " theta:" << options->getThetaobj() << endl;
      opt << " phi:" << options->getPhiobj() << endl;
      opt << " psi:" << options->getPsiobj() << endl;
    }
    else if (options->getObject() == "ellipsoid") {
      opt << " half axe A:" << options->getDemiaxea() << endl;
      opt << " half axe B:" << options->getDemiaxeb() << endl;
      opt << " half axe C:" << options->getDemiaxec() << endl;
      opt << " position X:" << options->getPositionx().at(i) << endl;
      opt << " position Y:" << options->getPositiony().at(i) << endl;
      opt << " position Z:" << options->getPositionz().at(i) << endl;
    }
    else if (options->getObject() == "cylinder") {
      opt << " radius:" << options->getSphereradius().at(i) << endl;
      opt << " height:" << options->getHauteur() << endl;
      opt << " position X:" << options->getPositionx().at(i) << endl;
      opt << " position Y:" << options->getPositiony().at(i) << endl;
      opt << " position Z:" << options->getPositionz().at(i) << endl;
      opt << " theta:" << options->getThetaobj() << endl;
      opt << " phi:" << options->getPhiobj() << endl;
    }
  }
    opt << " anisotropy:" << options->getAnisotropy() << endl;
  for (int i = 0 ; i < options->getObjectNumber(); i++) {
    if ( options->getAnisotropy() == "iso" ) {
      opt << "  material" << i << ":" << options->getMaterial()[i] << endl;
      opt << "  epsilon real:" << QString::number(real(options->getEpsilon().at(i))) << endl;
      opt << "  epsilon imag:" << QString::number(imag(options->getEpsilon().at(i))) << endl;
    }
    else if ( options->getAnisotropy() == "ani" ) {
      opt << "  epsilon11 real:" << QString::number(real(options->getEpsilon11().at(i))) << endl;
      opt << "  epsilon11 imag:" << QString::number(imag(options->getEpsilon11().at(i))) << endl;
      opt << "  epsilon12 real:" << QString::number(real(options->getEpsilon12().at(i))) << endl;
      opt << "  epsilon12 imag:" << QString::number(imag(options->getEpsilon12().at(i))) << endl;
      opt << "  epsilon13 real:" << QString::number(real(options->getEpsilon13().at(i))) << endl;
      opt << "  epsilon13 imag:" << QString::number(imag(options->getEpsilon13().at(i))) << endl;
      opt << "  epsilon21 real:" << QString::number(real(options->getEpsilon21().at(i))) << endl;
      opt << "  epsilon21 imag:" << QString::number(imag(options->getEpsilon21().at(i))) << endl;
      opt << "  epsilon22 real:" << QString::number(real(options->getEpsilon22().at(i))) << endl;
      opt << "  epsilon22 imag:" << QString::number(imag(options->getEpsilon22().at(i))) << endl;
      opt << "  epsilon23 real:" << QString::number(real(options->getEpsilon23().at(i))) << endl;
      opt << "  epsilon23 imag:" << QString::number(imag(options->getEpsilon23().at(i))) << endl;
      opt << "  epsilon31 real:" << QString::number(real(options->getEpsilon31().at(i))) << endl;
      opt << "  epsilon31 imag:" << QString::number(imag(options->getEpsilon31().at(i))) << endl;
      opt << "  epsilon32 real:" << QString::number(real(options->getEpsilon32().at(i))) << endl;
      opt << "  epsilon32 imag:" << QString::number(imag(options->getEpsilon32().at(i))) << endl;
      opt << "  epsilon33 real:" << QString::number(real(options->getEpsilon33().at(i))) << endl;
      opt << "  epsilon33 imag:" << QString::number(imag(options->getEpsilon33().at(i))) << endl;
    }
  }
  opt << "study options:" << endl;
  opt << " dipole/epsilon checked:" << options->getDipolepsilon() << endl;
  opt << " farfield checked:" << options->getFarfield() << endl;
  if ( options->getFarfield() ) {
    opt << "  cross section checked:" << options->getCrosssection() << endl;
    opt << "  cross section + poynting checked:" << options->getCrosssectionpoynting() << endl;
    opt << "  quick computation:" << options->getQuickdiffract() << endl;
    opt << "   ntheta:" << options->getNtheta() << endl;
    opt << "   nphi:" << options->getNphi() << endl;
    opt << "  microscopy checked:" << options->getMicroscopy() << endl;
    opt << "   numerical aperture:" << options->getNA() << endl;
    opt << "   magnification:" << options->getGross() << endl;
  }
  opt << " force checked:" << options->getForce() << endl;
  if ( options->getForce() ) {
    opt << "  optical force checked:" << options->getOpticalforce() << endl;
    opt << "  optical force density checked:" << options->getOpticalforcedensity() << endl;
    opt << "  optical torque checked:" << options->getOpticaltorque() << endl;
    opt << "  optical torque density checked:" << options->getOpticaltorquedensity() << endl;
  }
  opt << " nearfield checked:" << options->getNearfield() << endl;
  if ( options->getNearfield() ) {
    opt << "  localfield checked:" << options->getLocalfield() << endl;
    opt << "  macroscopic field checked:" << options->getMacroscopicfield() << endl;
    int nproche = options->getNproche() - 1;
    opt << "  range of study (-1,0,1):" << nproche << endl;
  }
  opt << "numerical parameters options:" << endl;
  opt << " tolerance:" << options->getTolerance() << endl;
  opt << " methode:" << options->getMethodeit() << endl;
  opt << " nxm:" << options->getNxm() << endl;
  opt << " nym:" << options->getNym() << endl;
  opt << " nzm:" << options->getNzm() << endl;
  opt << " polarizability:" << options->getPolarizability() << endl;
  opt << " quad:" << options->getQuad() << endl;
  opt << " FFT:" << options->getnfft2d() << endl;
  optfile.close();
  initOptionsTbl();
}
void 
OptionsWindow::showWarning(QString message) {
  QMessageBox::warning(this, "Error:", message);
}

