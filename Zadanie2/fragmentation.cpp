#include "fragmentation.h"

const double g_l1_max = 12;
const double g_l2_max = 12;
const double g_l1_min = 8;
const double g_l2_min = 8;
const double g_l0 = 5;

const double g_precision = 0.15;

/// вектор, содержащий box-ы, €вл€ющиес€ частью рабочего пространства
 std::vector<Box> solution;
/// вектор, содержащий box-ы, не €вл€ющиес€ частью рабочего пространства
std::vector<Box> not_solution;
/// вектор, содержащий box-ы, наход€щиес€ на границе между "рабочим" и "нерабочим" пространством
std::vector<Box> boundary;
/// вектор, хран€щий box-ы, анализируемые на следующей итерации алгоритма
std::vector<Box> temporary_boxes;

/// функции gj()
//------------------------------------------------------------------------------------------
double g1(double x1, double x2)
{
	return (x1*x1 + x2*x2 - g_l1_max*g_l1_max);
}

//------------------------------------------------------------------------------------------
double g2(double x1, double x2)
{
	return (g_l1_min*g_l1_min - x1*x1 - x2*x2);
}

//------------------------------------------------------------------------------------------
double g3(double x1, double x2)
{
	return (x1*x1 + x2*x2 - g_l2_max*g_l2_max);
}

//------------------------------------------------------------------------------------------
double g4(double x1, double x2)
{
	return (g_l2_min*g_l2_min - x1*x1 - x2*x2);
}


//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(double& min_x, double& min_y, double& x_width, double& y_height )
{
	current_box = Box( min_x, min_y, x_width, y_height );
}

//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(const Box& box)
{
	current_box = box;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::VerticalSplitter(const Box& box, boxes_pair& vertical_splitter_pair)
{
	// необходимо определить функцию.
	double x, y, w, h;
	box.GetParameters(x, y, w, h);
	Box b1 = Box(x, y, w/2, h);
	Box b2 = Box(x + w/2, y, w / 2, h);

	vertical_splitter_pair.first = b1;
	vertical_splitter_pair.second = b2;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::HorizontalSplitter(const Box& box, boxes_pair& horizontal_splitter_pair)
{
	// необходимо определить функцию
	double x, y, w, h;
	box.GetParameters(x, y, w, h);
	Box b1 = Box(x, y, w, h/2);
	Box b2 = Box(x, y + h/2, w, h/2);

	horizontal_splitter_pair.first = b1;
	horizontal_splitter_pair.second = b2;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes)
{
	// необходимо определить функцию
	double x, y, w, h;
	box.GetParameters(x, y, w, h);
	if (w > h) {
		VerticalSplitter(box, new_pair_of_boxes);
	}
	else
	{
		HorizontalSplitter(box, new_pair_of_boxes);
	}
}

//------------------------------------------------------------------------------------------
unsigned int low_level_fragmentation::FindTreeDepth()
{
	double box_diagonal = current_box.GetDiagonal();

	if (box_diagonal <= g_precision)
	{
		return 0;
	}
	else
	{
		boxes_pair new_boxes;
		// допустим, разобьем начальную область по ширине
		VerticalSplitter(current_box, new_boxes);
		unsigned int tree_depth = 1;

		box_diagonal = new_boxes.first.GetDiagonal();

		if (box_diagonal <= g_precision)
		{
			return tree_depth;
		}
		else
		{
			for (;;)
			{
				GetNewBoxes(new_boxes.first, new_boxes);
				++tree_depth;
				box_diagonal = new_boxes.first.GetDiagonal();

				if (box_diagonal <= g_precision)
				{
					break;
				}
			}
			return tree_depth;
		}
	}
}

//------------------------------------------------------------------------------------------
int low_level_fragmentation::ClasifyBox(const min_max_vectors& vects)
{
	// необходимо определить функцию

	bool max_cond = true;
	std::cout << "max" << std::endl;
	for (double el:vects.second) {
		std::cout << "class " << el << std::endl;
		if (el>0) {
			max_cond = false;
			//break;
		}
	}

	//contain points
	if (max_cond) return 0;

	std::cout << "min" << std::endl;
	bool min_cond = false;
	int c = 0;
	for (double el : vects.first) {
		if (el > 0) {
			std::cout << "class " << el << std::endl;
			c++;
			//min_cond = true;
			break;
		}
	}
	if (c > 0) min_cond = true;
	if (min_cond) return 1;

	//nor
	return 2;
}

//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetBoxType(const Box& box)
{
	// необходимо определить функцию
	min_max_vectors verts;
	GetMinMax(box, verts);
	int clas = ClasifyBox(verts);
	//std::cout << "class " << clas << std::endl;
	switch (clas)
	{
	case 0:
		//std::cout << "class " << clas << std::endl;
		boundary.push_back(box);
		break;

	case 1:
		std::cout << "class " << clas << std::endl;
		not_solution.push_back(box);
		break;

	case 2:
		boxes_pair boxes;
		solution.push_back(box);
		GetNewBoxes(box, boxes);
		temporary_boxes.push_back(boxes.first);
		temporary_boxes.push_back(boxes.second);
		break;
	}
}


//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis( double& min_x, double& min_y, double& x_width, double& y_height ) :
					low_level_fragmentation(min_x, min_y, x_width, y_height) {}

//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis( Box& box ) : low_level_fragmentation( box ) {}

//------------------------------------------------------------------------------------------
void high_level_analysis::GetMinMax( const Box& box, min_max_vectors& min_max_vecs )
{
	std::vector<double> g_min;
	std::vector<double> g_max;

	double a1min, a2min, a1max, a2max;
	double xmin, xmax, ymin, ymax;

	box.GetParameters(xmin, ymin, xmax, ymax);

	xmax = xmin + xmax;
	ymax = ymin + ymax;

	double curr_box_diagonal = box.GetDiagonal();

	if (curr_box_diagonal <= g_precision)
	{
		g_min.push_back(0);
		g_max.push_back(0);

		min_max_vecs.first = g_min;
		min_max_vecs.second = g_max;

		return;
	}

	// MIN
	// функци€ g1(x1,x2)
	a1min = __min(abs(xmin), abs(xmax));
	a2min = __min(abs(ymin), abs(ymax));
	g_min.push_back(g1(a1min, a2min));

	// функци€ g2(x1,x2)
	a1min = __max(abs(xmin), abs(xmax));
	a2min = __max(abs(ymin), abs(ymax));
	g_min.push_back(g2(a1min, a2min));

	// функци€ g3(x1,x2)
	a1min = __min(abs(xmin - g_l0), abs(xmax - g_l0));
	a2min = __min(abs(ymin), abs(ymax));
	g_min.push_back(g3(a1min, a2min));

	// функци€ g4(x1,x2)
	a1min = __max(abs(xmin - g_l0), abs(xmax - g_l0));
	a2min = __max(abs(ymin), abs(ymax));
	g_min.push_back(g4(a1min, a2min));

	// MAX
	// функци€ g1(x1,x2)
	a1max = __max(abs(xmin), abs(xmax));
	a2max = __max(abs(ymin), abs(ymax));
	g_max.push_back(g1(a1max, a2max));

	// функци€ g2(x1,x2)
	a1max = __min(abs(xmin), abs(xmax));
	a2max = __min(abs(ymin), abs(ymax));
	g_max.push_back(g2(a1max, a2max));

	// функци€ g3(x1,x2)
	a1max = __max(abs(xmin - g_l0), abs(xmax - g_l0));
	a2max = __max(abs(ymin), abs(ymax));
	g_max.push_back(g3(a1max, a2max));

	// функци€ g4(x1,x2)
	a1max = __min(abs(xmin - g_l0), abs(xmax - g_l0));
	a2max = __min(abs(ymin), abs(ymax));
	g_max.push_back(g4(a1max, a2max));

	min_max_vecs.first = g_min;
	min_max_vecs.second = g_max;
}

//------------------------------------------------------------------------------------------
void high_level_analysis::GetSolution()
{
	// необходимо определить функцию
	double cur_sigma = current_box.GetDiagonal();

	int tree_depth = FindTreeDepth();
	//std::cout << "kek " <<  tree_depth <<  std::endl;
	
	solution.push_back(current_box);
	temporary_boxes.push_back(current_box);

	for (int i = 0; i < tree_depth;++i) {
		std::vector<Box> tmp = temporary_boxes;
		//std::cout << "Num " << tmp.size() << std::endl;
		temporary_boxes.clear();
		for(Box b:tmp)
		{
			GetBoxType(b);
		}
	}
}


//------------------------------------------------------------------------------------------
void WriteResults( const char* file_names[] )
{
	// необходимо определить функцию

	//in
	std::ofstream myfile;
	myfile.open(file_names[0]);
	for (Box b:solution) {
		//std::cout << "kek" << std::endl;
		double x, y, w, h;
		b.GetParameters(x, y, w, h);
		myfile << x << ", " << y << ", " << w << ", " << h << ", \n";
	}
	myfile.close();

	//out
	myfile.open(file_names[1]);
	for (Box b : not_solution) {
		//std::cout << "kek" << std::endl;
		double x, y, w, h;
		b.GetParameters(x, y, w, h);
		myfile << x << ", " << y << ", " << w << ", " << h << ", \n";
	}
	myfile.close();

	//bound
	myfile.open(file_names[2]);
	for (Box b : boundary) {
		//std::cout << "kek" << std::endl;
		double x, y, w, h;
		b.GetParameters(x, y, w, h);
		myfile << x << ", " << y << ", " << w << ", " << h << ", \n";
	}
	myfile.close();
}