

double toDoubleAdv(char * argv[], int argc, int & idx, char const * name)
{
	if (idx >= argc)
	{
		throw std::exception(QString("Expected more arguments, only %1 given, param %2\n").arg(argc).arg(name).toStdString().c_str());
	}
	bool ok = true;
	double result = QString(argv[idx++]).toDouble(&ok);
	if (!ok)
	{
		throw std::exception(QString("Argument parsing: Conversion of %1 to double failed, param %2!\n").arg(argv[idx - 1]).arg(name).toStdString().c_str());
	}
	return result;
}
int toIntAdv(char * argv[], int argc, int & idx, char const * name)
{
	if (idx >= argc)
	{
		throw std::exception(QString("Expected more arguments, only %1 given, param %2\n").arg(argc).arg(name).toStdString().c_str());
	}
	bool ok = true;
	int result = QString(argv[idx++]).toInt(&ok);
	if (!ok)
	{
		throw std::exception(QString("Argument parsing: Conversion of %1 to int failed, param %2!\n").arg(argv[idx - 1]).arg(name).toStdString().c_str());
	}
	return result;
}

bool isDouble(QString const & str)
{
	bool ok;
	str.toDouble(&ok);
	return ok;
}