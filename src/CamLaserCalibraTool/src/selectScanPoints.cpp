//
// Created by heyijia on 18-12-11.
//
#include "selectScanPoints.h"

// int img_w = 720;
// double focal = 450;
// double z = 10;    // 在 10m 高处装一个相机，咔咔给激光点云拍照

int img_w = 1280;
double focal = 15;
double z = 0.2;    // 在 0.2m 高处装一个相机，咔咔给激光点云拍照

static std::mutex g_csv_mtx;
static std::ofstream g_csv;
static bool g_csv_inited = false;

struct LineSeg
{
    int id_start;
    int id_end;
    double dist;
};

std::vector< Eigen::Vector3d > AutoGetLinePts(const std::vector<Eigen::Vector3d> points, bool debug)
{
//    cv::Mat img(img_w, img_w, CV_8UC1, cv::Scalar::all(0));
    std::lock_guard<std::mutex> lk(g_csv_mtx);
    if (!g_csv_inited) {
        g_csv.open("/home/zhu/LaserCameraCal_ws/data/laser_points_debug.csv", std::ios::out | std::ios::app);
        if (g_csv.is_open()) {
            g_csv << "dist,len,count\n";
            g_csv_inited = true;
             printf("open csv!\n");
        }
        printf("init csv\n!");
    }

    cv::Mat img(img_w, img_w, CV_8UC3, cv::Scalar(0,0,0));

    for (auto pt: points) {
        int col = (int)(pt.x() / z * focal + img_w/2);
        int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果

        if(col > img_w-1 || col< 0 || row > img_w-1 || row < 0)
            continue;

        cv::Vec3b color_value(255,0,0);
        img.at<cv::Vec3b>(row, col) = color_value;
//        img.at<uchar>(row, col) = 255;
    }

    /// detect and get line
    // 直接从每一帧激光的正前方开始搜索一定距离范围内的符合平面标定板形状的激光线段
    int n = static_cast<int>(points.size());
    if (n==0) return {};
    int id = -1000;
//        std::cout << points.at(id).transpose() <<" "<<points.at(id+1).transpose() <<std::endl;
    // 假设每个激光点之间的夹角为0.3deg,
    // step 1: 如果有激光标定板，那么激光标定板必须出现在视野的正前方 120 deg 范围内(通常相机视野也只有 120 deg)，也就是左右各 60deg.
    int delta = 120/0.1;
    // const int delta = static_cast<int>(std::round(80.0/0.3));

    // int id_left = std::min( id + delta, n-1);  // 200
    // int id_right = std::max( id - delta , 0);  // 0
    int id_left = n-1;
    int id_right= 0;

    double dist_left  = points.at(id_left).head<2>().norm();
    double dist_right = points.at(id_right).head<2>().norm();
    // printf("AutoGetLinePts n: %d\n", n);
    // printf("AutoGetLinePts id_left: %d, id_right: %d\n", id_left, id_right);
    // printf("AutoGetLinePts dist_left: %.3f, dist_right: %.3f\n", dist_left, dist_right);

    // 边界检查
    auto inRange = [&](int k){ return k >= 0 && k < n; };
    auto range2  = [&](int k){ const auto& p = points[k]; return p.x()*p.x()+p.y()*p.y(); };

    // 逻辑别搞复杂了。
    std::vector<LineSeg> segs;
    double dist_thre = 0.004;
    int skip = 1;
    const double range_max2 = 100.0*100.0;
    int currentPt = id_right; // 0
    int nextPt = currentPt + skip; // 1
    bool newSeg = true;
    LineSeg seg;
    // for (int i = id_right; i < id_left - skip; i += skip) {
    for (int i = id_right; i <= id_left; i += skip) {
        if (!inRange(currentPt)) break;
        nextPt = currentPt + skip;
        if (!inRange(nextPt)) break;

        if(newSeg)
        {
            seg.id_start = currentPt;
            seg.id_end = nextPt;
            newSeg = false;
        }

        // head(2) 取坐标点前两个值x、y, norm() 取欧式距离
        // double d1 = points.at(currentPt).head(2).norm();
        // double d2 = points.at(nextPt).head(2).norm();

        double d1_2 = range2(currentPt);
        double d2_2 = range2(nextPt);
        if (!std::isfinite(d1_2) || !std::isfinite(d2_2)) { currentPt = nextPt; continue; }

        // double range_max = 100;
        if (d1_2 < range_max2 && d2_2 < range_max2) {
            // 与原逻辑一致：用半径差阈值
            // 用相邻点的“半径差”判断是否属于同一线段：对于平面靶反射的点，它们到雷达的距离在短邻域内变化不大
            double d1 = std::sqrt(d1_2), d2 = std::sqrt(d2_2);
            // printf("d1-d2=%lf\n", d1-d2);
            if (std::fabs(d1 - d2) < dist_thre) {
                seg.id_end = nextPt;
            } else {
                newSeg = true;
                if (inRange(seg.id_start) && inRange(seg.id_end)) {
                    Eigen::Vector2d dd = (points[seg.id_start] - points[seg.id_end]).head<2>();
                    double len = dd.norm();  // 求解该段的几何长度（平面
                    // printf("len=%lf\n", len);
                    // printf("seg.id_end - seg.id_start=%d\n", seg.id_end - seg.id_start);
                    if (g_csv.is_open()) {
                        const double dist = d1 - d2;
                        const int count = seg.id_end - seg.id_start;
                        g_csv << dist << ',' << len << ',' << count << '\n';
                        g_csv.flush();
                    }
                    
                    if ( len > 0.02
                        && std::sqrt(range2(seg.id_start)) < 2  // 两端点都在2 m内（更可能是近处靶）
                        && std::sqrt(range2(seg.id_end))   < 2
                        && (seg.id_end - seg.id_start) > 100)  // 点数>50（稠密，稳定）
                    {
                        seg.dist = len;
                        segs.push_back(seg);
                    }
                }
            }
            currentPt = nextPt;
        } else {
            currentPt = nextPt;
        }
    }


    // 对 segs 的边界进行扩充
    // for (int i = 0; i < segs.size(); ++i) {
    //     LineSeg tmp = segs.at(i);

    //     for (int j = 1; j < 4; ++j)
    //     {
    //         int boundaryPt = tmp.id_end + j;
    //         double d1 = points.at(tmp.id_end).head(2).norm();
    //         double d2 = points.at(boundaryPt).head(2).norm();
    //         if(fabs(d1-d2) < dist_thre)  //  8cm
    //         {
    //             segs.at(i).id_end = boundaryPt;
    //         }
    //     }

    //     for (int j = -1; j > -4; --j)
    //     {
    //         int boundaryPt = tmp.id_start + j;
    //         double d1 = points.at(tmp.id_start).head(2).norm();
    //         double d2 = points.at(boundaryPt).head(2).norm();
    //         if(fabs(d1-d2) < dist_thre)  //  8cm
    //         {
    //             segs.at(i).id_start = boundaryPt;
    //         }
    //     }

    // }

    // ——边界扩展：加上越界保护；一旦越界或不满足条件就停止扩展
    for (auto& s : segs) {

        // 右扩：+1..+3
        for (int j = 1; j <= 3; ++j) {
            int b = s.id_end + j;
            if (!inRange(b)) break;
            double r_end = std::sqrt(range2(s.id_end));
            double r_b   = std::sqrt(range2(b));
            if (!std::isfinite(r_end) || !std::isfinite(r_b)) break;
            if (std::fabs(r_end - r_b) < dist_thre) s.id_end = b; else break;
        }

        // 左扩：-1..-3（关键修复点）
        for (int j = 1; j <= 3; ++j) {
            int b = s.id_start - j;
            if (!inRange(b)) break;                          // ✅ 防止 -1
            double r_sta = std::sqrt(range2(s.id_start));
            double r_b   = std::sqrt(range2(b));
            if (!std::isfinite(r_sta) || !std::isfinite(r_b)) break;
            if (std::fabs(r_sta - r_b) < dist_thre) s.id_start = b; else break;
        }
    }

    // std::vector<Eigen::Vector3d> ptsLine;
    // if(segs.size() > 0)
    // {
    //     LineSeg bestLine;
    //     int maxpts = -1;
    //     for (int i = 0; i < segs.size(); ++i) {
    //         LineSeg tmp = segs.at(i);
    //         int cnt = tmp.id_end - tmp.id_start;
    //         if(cnt > maxpts)
    //         {
    //             bestLine = tmp;
    //             maxpts = cnt;
    //         }
    //     }

    //     for (int i = bestLine.id_start; i < bestLine.id_end+1; ++i) {
    //         ptsLine.push_back(points.at(i));
    //     }
    // }

    // ——选最长线段：同样做范围检查
    std::vector<Eigen::Vector3d> ptsLine;
    if (!segs.empty()) {
        auto best = std::max_element(segs.begin(), segs.end(),
                                    [](const LineSeg& a, const LineSeg& b){
                                        return (a.id_end - a.id_start) < (b.id_end - b.id_start);
                                    });
        for (int i = best->id_start; i <= best->id_end; ++i)
            if (inRange(i)) ptsLine.push_back(points[i]);
    }


    printf("debug: %d\n", debug);
    if(debug)
    {
        for (int j = 0; j < segs.size(); ++j) {
            LineSeg tmp = segs.at(j);
            for (int i = tmp.id_start; i < tmp.id_end; ++i)
            {
                Eigen::Vector3d pt = points.at(i);
                int col = (int)(pt.x() / z * focal + img_w/2);
                int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果

                if(col > img_w-1 || col< 0 || row > img_w-1 || row < 0)
                    continue;

                cv::Vec3b color_value(0,255,0);
                img.at<cv::Vec3b>(row, col) = color_value;
            }
        }

        for (int j = 0; j < ptsLine.size(); ++j) {

            Eigen::Vector3d pt = ptsLine.at(j);
            int col = (int)(pt.x() / z * focal + img_w/2);
            int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果

            if(col > img_w-1 || col< 0 || row > img_w-1 || row < 0)
                continue;

            cv::Vec3b color_value(0, 0, 255); // BGR顺序
            img.at<cv::Vec3b>(row, col) = color_value;

        }

        cv::putText(img, "Detecting the Laser Points on the calibra planar!",
                    cv::Point(5,30), cv::FONT_HERSHEY_COMPLEX_SMALL, 0.7, cv::Scalar(255,255,255), 1, cv::LINE_AA);
        cv::imshow("ScanPoint",img);
        cv::waitKey(10);
    }

    return ptsLine;

}

std::vector<cv::Rect> PointsToImg(const std::vector<Eigen::Vector3d> points, bool selectPoint)
{

    cv::Mat img(img_w, img_w, CV_8UC1, cv::Scalar::all(0));

    for (auto pt: points) {
        int col = (int)(pt.x() / z * focal + img_w/2);
        int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果

        if(col > img_w-1 || col< 0 || row > img_w-1 || row < 0)
            continue;

        img.at<uchar>(row, col) = 255;
    }

    std::vector<cv::Rect> rects;
    if(selectPoint)
    {
        GetRect getrect;
        getrect.gettingROI(img);  // 选取感兴趣的图像区域
        rects = getrect.rects;
    } else
    {
        cv::imshow("ScanPoint",img);
        cv::waitKey(0);
    }

    return rects;

}


std::vector<Eigen::Vector3d> GetROIScanPoints(const std::vector<Eigen::Vector3d> points, const std::vector<cv::Rect> rects)
{
    std::vector<Eigen::Vector3d> RoiPoints;
    for (auto rect : rects) {

        for (auto pt: points)
        {
            // 将点转换成图像坐标，然后看是否落在感兴趣区域内部
            int col = (int)(pt.x() / z * focal + img_w/2);
            int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果
            if(col > rect.x + rect.width || col< rect.x || row > rect.y +rect.height || row < rect.y)
                continue;
            RoiPoints.push_back(pt);
        }
    }
    return RoiPoints;
}

void GetRect::onMouse(int event, int x, int y, int flags, void *userdata) {
    // Check for null pointer in userdata and handle the error
    GetRect* temp = reinterpret_cast<GetRect*>(userdata);
    temp->CallBackFunc(event, x, y);

}
void GetRect::CallBackFunc(int event, int x, int y) {

    if ( event == cv::EVENT_LBUTTONDOWN )
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;

        // Init your rect
        base.x = x;
        base.y = y;
        r.x = x;
        r.y = y;
        r.width = 0;
        r.height = 0;
        bDraw = true;
    }
    else if ( event == cv::EVENT_MOUSEMOVE )
    {
//        std::cout << "Mouse move over the window - position (" << x << ", " << y << ")" << std::endl;

        // If drawing, update rect width and height
        if(!bDraw) return;

        int dx = abs(r.x - x);
        int dy = abs(r.y - y);

        if(x < base.x) {
            r.x = x;
            r.width = abs(x - base.x);
        } else {
            r.width = dx;
        }

        if(y < base.y) {
            r.y = y;
            r.height = abs(y - base.y);
        } else {
            r.height = dy;
        }

        // Refresh
        working = layer.clone();
        cv::rectangle(working, r, cv::Scalar(125));
        cv::imshow("ScanPoint", working);
    }
    else if ( event == cv::EVENT_LBUTTONUP)
    {
        std::cout << "Left button released" << std::endl;

        // Save rect, draw it on layer
        rects.push_back(r);
        cv::rectangle(layer, r, cv::Scalar(125));

        r = cv::Rect();
        bDraw = false;

        // Refresh
        working = layer.clone();
        cv::rectangle(working, r, cv::Scalar(125));
        cv::imshow("ScanPoint", working);
    }
}
void GetRect::gettingROI(cv::Mat img)
{
    layer = img.clone();
    cv::putText(layer, "Please Select ROI scan data with your mouse. Then, Put any key to continue!",
                cv::Point(5,30), cv::FONT_HERSHEY_COMPLEX_SMALL, 0.7, cv::Scalar(255), 1, cv::LINE_AA);
    cv::namedWindow("ScanPoint", 1);
    cv::imshow("ScanPoint", layer);
    cv::setMouseCallback("ScanPoint", GetRect::onMouse, this);
    cv::waitKey(0);
}
