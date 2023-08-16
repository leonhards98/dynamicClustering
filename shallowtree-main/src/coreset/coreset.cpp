#include "coreset.hpp"
#include <map>

using namespace coresets;

Coreset::Coreset(size_t targetSize, utils::IndexFinder &i) : TargetSize(targetSize), i(i)
{
}


void Coreset::writePointsToVector(std::vector<dynClustering::CoresetPoint*>* pointVector, std::vector<dynClustering::CoresetPoint*>* centersVector)
{
    pointVector->clear();
    centersVector->clear();
    int counter {-1};
    for (auto p = this->points.begin(); p != this->points.end(); p++)
    {
        if ((*p)->IsCenter)
        {
            counter++;
			dynClustering::CoresetPoint* c = new dynClustering::CoresetPoint(-1,(*p)->Weight,*this->centers.at(counter));
            centersVector->push_back(c);
            pointVector->push_back(c);
            
            continue;
        }

        pointVector->push_back(this->i.getPoint((*p)->Index));
    }
}

void Coreset::writeWeightsToVector(std::vector<double>* weightVector)
{
    weightVector->clear();
    for (auto p = this->points.begin(); p != this->points.end(); p++)
    {
        weightVector->push_back((*p)->Weight);
    }
}

std::shared_ptr<WeightedPoint>
Coreset::findPoint(size_t index, bool isCenter)
{
    for (size_t i = 0; i < this->points.size(); i++)
    {
        if (points[i]->Index == index && points[i]->IsCenter == isCenter)
        {
            return points[i];
        }
    }

    return nullptr;
}

void Coreset::addPoint(size_t pointIndex, double weight)
{
    //auto coresetPoint = findPoint(pointIndex, false);
    //if (coresetPoint == nullptr)
    //{
        //coresetPoint = std::make_shared<WeightedPoint>(pointIndex, 0.0, false);
        auto coresetPoint = std::make_shared<WeightedPoint>(pointIndex, weight, false);
        this->points.push_back(coresetPoint);
    //}
    //else
    //{
        // Updating
    //}

    //coresetPoint->Weight += weight;
}

void Coreset::addCenter(size_t clusterIndex, std::shared_ptr<blaze::DynamicVector<double>> center, double weight)
{
    auto coresetPoint = findPoint(clusterIndex, true);
    if (coresetPoint == nullptr)
    {
        coresetPoint = std::make_shared<WeightedPoint>(clusterIndex, 0.0, true);
        this->points.push_back(coresetPoint);
        centers.emplace(clusterIndex, center);
    }
    else
    {
        // Updating
    }
    coresetPoint->Weight += weight;
}

std::shared_ptr<WeightedPoint>
Coreset::at(size_t index) const
{
    return this->points[index];
}

size_t
Coreset::size() const
{
    return this->points.size();
}
void Coreset::writeToStream(std::ostream &out)
{
    out << "Points:\n";
    for (auto x : this->points)
    {
        out << x->Index << "\n";
        
    }

    out << "Centers:\n";
    for (auto x : this->centers)
    {
        out << *x.second << "\n";
    }
}

void Coreset::writeCoordsToStream(std::ostream &out)
{
    for (auto x : this->points)
    {
        if (x->IsCenter)
        {
            out << "Center: " << x->Weight << " " << *this->centers.at(x->Index) << "\n";
            continue;
        }
        out << "Point: " << x->Weight << " " << this->i.getPoint(x->Index)->getCoordinates() << "\n";
        
    }

}

void Coreset::writeToStream(const blaze::DynamicMatrix<double> &originalDataPoints, std::ostream &out)
{
    const size_t m = this->points.size();
    const size_t d = originalDataPoints.columns();
	size_t centers {0};

    // Output coreset size
    out << m << "\n";

    std::shared_ptr<blaze::DynamicVector<double>> center;

    // Output coreset points
    for (auto &&point : points)
    {
        // Output coreset point weight
        out << point->Weight << " ";

        if (point->IsCenter)
        {
			out << "Center: ";
            center = this->centers.at(point->Index);
			centers++;
        }

        // Output coreset point entries.
        for (size_t j = 0; j < d; ++j)
        {
            if (point->IsCenter)
            {
                out << center->at(j);
            }
            else
            {
                out << originalDataPoints.at(point->Index, j);
            }

            if (j < d - 1)
            {
                out << " ";
            }
        }
        out << "\n";
    }
	out << "Centers: " << centers << "\n";
}
